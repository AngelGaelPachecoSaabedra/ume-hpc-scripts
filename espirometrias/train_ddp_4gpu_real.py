# -*- coding: utf-8 -*-
import os, re, json, argparse, random
from typing import List, Tuple, Optional
import numpy as np
import pandas as pd

import torch
import torch.nn as nn
import torch.distributed as dist
from torch.nn.parallel import DistributedDataParallel as DDP
from torch.utils.data import Dataset, DataLoader
from torch.utils.data.distributed import DistributedSampler
from torchvision import transforms, models
from PIL import Image, UnidentifiedImageError
from sklearn.model_selection import GroupShuffleSplit
from tqdm import tqdm

# ==================== DDP Setup ====================

def setup_ddp():
    """Inicializar proceso distribuido"""
    dist.init_process_group(backend='nccl')
    rank = int(os.environ.get('RANK', 0))
    local_rank = int(os.environ.get('LOCAL_RANK', 0))
    world_size = int(os.environ.get('WORLD_SIZE', 1))
    torch.cuda.set_device(local_rank)
    return rank, local_rank, world_size

# ==================== Utilidades ====================

def seed_everything(seed=42):
    random.seed(seed); np.random.seed(seed)
    torch.manual_seed(seed); torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

def split_paths(cell: str) -> List[str]:
    if pd.isna(cell): return []
    return [p.strip() for p in str(cell).split("|") if str(p).strip() != ""]

def product_pairs(front_list: List[str], side_list: List[str]) -> List[Tuple[str,str]]:
    if not front_list or not side_list: return []
    return [(f, s) for f in front_list for s in side_list]

def load_rgb(path: str) -> Image.Image:
    return Image.open(path).convert("RGB")

def _to_num(x):
    if pd.isna(x): return np.nan
    s = str(x).strip().lower()
    s = s.replace(",", ".")
    s = re.sub(r"[^0-9.\-eE+]", "", s)
    if s in {"", ".", "-", "+", "e", "E"}: return np.nan
    try:
        return float(s)
    except Exception:
        return np.nan

# ==================== Transforms ====================

class PairedAugment:
    def __init__(self, img_size=224):
        self.img_size = img_size
    def __call__(self, img_front: Image.Image, img_side: Image.Image):
        img_front = img_front.resize((self.img_size, self.img_size))
        img_side  = img_side.resize((self.img_size, self.img_size))
        if random.random() < 0.5:
            img_front = img_front.transpose(Image.FLIP_LEFT_RIGHT)
            img_side  = img_side.transpose(Image.FLIP_LEFT_RIGHT)
        if random.random() < 0.3:
            b=c=s=0.1; h=0.05
            def jitter(img):
                import torchvision.transforms.functional as F
                img = F.adjust_brightness(img, 1 + random.uniform(-b,b))
                img = F.adjust_contrast(img,  1 + random.uniform(-c,c))
                img = F.adjust_saturation(img, 1 + random.uniform(-s,s))
                img = F.adjust_hue(img,       random.uniform(-h,h))
                return img
            img_front = jitter(img_front); img_side = jitter(img_side)
        to_tensor = transforms.ToTensor()
        norm = transforms.Normalize(mean=[0.485,0.456,0.406], std=[0.229,0.224,0.225])
        return norm(to_tensor(img_front)), norm(to_tensor(img_side))

class EvalTransform:
    def __init__(self, img_size=224):
        self.tf = transforms.Compose([
            transforms.Resize((img_size,img_size)),
            transforms.ToTensor(),
            transforms.Normalize(mean=[0.485,0.456,0.406], std=[0.229,0.224,0.225]),
        ])
    def __call__(self, img_front, img_side):
        return self.tf(img_front), self.tf(img_side)

# ==================== Dataset ====================

class DuoViewsImgTabDataset(Dataset):
    def __init__(self, df, tab_cols: List[str], train=True, img_size=224,
                 y_mean=None, y_std=None, tab_mean=None, tab_std=None):
        self.df = df.reset_index(drop=True)
        self.train = train
        self.aug = PairedAugment(img_size) if train else EvalTransform(img_size)
        self.tab_cols = tab_cols

        y = self.df[["peso_kg","altura_cm","ict"]].values.astype(np.float32)
        if y_mean is None or y_std is None:
            self.y_mean = y.mean(axis=0, keepdims=True)
            self.y_std  = y.std(axis=0, keepdims=True) + 1e-6
        else:
            self.y_mean, self.y_std = y_mean, y_std
        self.y_norm = (y - self.y_mean) / self.y_std

        if len(self.tab_cols):
            Xtab = self.df[self.tab_cols].values.astype(np.float32)
            if tab_mean is None or tab_std is None:
                self.tab_mean = Xtab.mean(axis=0, keepdims=True)
                self.tab_std  = Xtab.std(axis=0, keepdims=True) + 1e-6
            else:
                self.tab_mean, self.tab_std = tab_mean, tab_std
            self.tab_norm = (Xtab - self.tab_mean) / self.tab_std
        else:
            self.tab_mean = None; self.tab_std = None
            self.tab_norm = None

    def __len__(self): return len(self.df)

    def __getitem__(self, idx):
        r = self.df.iloc[idx]
        fpath, spath = r["front_path"], r["side_path"]
        try:
            img_f = load_rgb(fpath); img_s = load_rgb(spath)
        except (FileNotFoundError, UnidentifiedImageError) as e:
            raise RuntimeError(f"No pude abrir:\nfront={fpath}\nside={spath}\n{e}")
        x_f, x_s = self.aug(img_f, img_s)
        y = torch.tensor(self.y_norm[idx], dtype=torch.float32)
        if self.tab_norm is None:
            x_tab = torch.empty(0, dtype=torch.float32)
        else:
            x_tab = torch.tensor(self.tab_norm[idx], dtype=torch.float32)
        return {"front": x_f, "side": x_s, "tab": x_tab, "y": y, "id": r["encuesta_id"]}

# ==================== Modelo ====================

class TwoBranchImgTabRegressor(nn.Module):
    def __init__(self, backbone="resnet18", dropout=0.3, out_dim=3, tab_dim: int = 0, tab_hidden=64):
        super().__init__()
        if backbone == "efficientnet_b0":
            from torchvision.models import EfficientNet_B0_Weights
            try:
                base = models.efficientnet_b0(weights=EfficientNet_B0_Weights.DEFAULT)
            except Exception:
                base = models.efficientnet_b0(pretrained=True)
            feat_dim = base.classifier[1].in_features
            base.classifier = nn.Identity()
            self.encoder = base
        elif backbone == "resnet18":
            from torchvision.models import ResNet18_Weights
            try:
                base = models.resnet18(weights=ResNet18_Weights.DEFAULT)
            except Exception:
                base = models.resnet18(pretrained=True)
            feat_dim = base.fc.in_features
            base.fc = nn.Identity()
            self.encoder = base
        elif backbone == "efficientnet_b2":
            from torchvision.models import EfficientNet_B2_Weights
            try:
                base = models.efficientnet_b2(weights=EfficientNet_B2_Weights.DEFAULT)
            except Exception:
                base = models.efficientnet_b2(pretrained=True)
            feat_dim = base.classifier[1].in_features
            base.classifier = nn.Identity()
            self.encoder = base
        elif backbone == "efficientnet_b3":
            from torchvision.models import EfficientNet_B3_Weights
            try:
                base = models.efficientnet_b3(weights=EfficientNet_B3_Weights.DEFAULT)
            except Exception:
                base = models.efficientnet_b3(pretrained=True)
            feat_dim = base.classifier[1].in_features
            base.classifier = nn.Identity()
            self.encoder = base
        else:
            raise ValueError("Backbone no soportado")

        self.tab_dim = int(tab_dim)
        if self.tab_dim > 0:
            self.tab_mlp = nn.Sequential(
                nn.Linear(self.tab_dim, tab_hidden),
                nn.ReLU(inplace=False),
                nn.Dropout(dropout),
                nn.Linear(tab_hidden, tab_hidden),
                nn.ReLU(inplace=False),
            )
            tab_out = tab_hidden
        else:
            self.tab_mlp = None
            tab_out = 0

        fusion_in = feat_dim*2 + tab_out
        self.head = nn.Sequential(
            nn.Linear(fusion_in, 512),
            nn.ReLU(inplace=False),
            nn.Dropout(dropout),
            nn.Linear(512, out_dim),
        )

    def forward(self, front, side, tab: Optional[torch.Tensor] = None):
        f = self.encoder(front)
        s = self.encoder(side)
        z = torch.cat([f, s], dim=1)
        if self.tab_mlp is not None and tab is not None and tab.nelement() > 0:
            t = self.tab_mlp(tab)
            z = torch.cat([z, t], dim=1)
        return self.head(z)

# ==================== Expansión vistas ====================

def expand_all_photos(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in df.iterrows():
        fronts = split_paths(r["fotos_frontal"])
        sides  = split_paths(r["fotos_lateral"])
        for (f,s) in product_pairs(fronts, sides):
            rows.append({
                "encuesta_id": r["encuesta_id"],
                "front_path": f,
                "side_path": s,
                "peso_kg": r["peso_kg"],
                "altura_cm": r["altura_cm"],
                "cintura_cm": r["cintura_cm"],
            })
    out = pd.DataFrame(rows)
    if out.empty: return out
    out = out[
        out["front_path"].apply(lambda p: os.path.exists(str(p))) &
        out["side_path"].apply(lambda p: os.path.exists(str(p)))
    ].reset_index(drop=True)
    out["ict"] = out["cintura_cm"] / np.clip(out["altura_cm"], 1e-6, None)
    return out

# ==================== Métricas ====================

def mae(a,b): return float(np.mean(np.abs(a-b)))
def rmse(a,b): return float(np.sqrt(np.mean((a-b)**2)))

# ==================== Entrenamiento ====================

def main(args):
    # Setup DDP
    rank, local_rank, world_size = setup_ddp()
    device = torch.device(f"cuda:{local_rank}")

    # Solo rank 0 imprime
    def print_rank0(*pargs, **kwargs):
        if rank == 0:
            print(*pargs, **kwargs, flush=True)

    seed_everything(args.seed + rank)  # Seed diferente por rank

    if rank == 0:
        os.makedirs(args.out_dir, exist_ok=True)

    # Leer CSV
    base = pd.read_csv(args.csv_path)
    rename_map = {}
    if "peso" in base.columns: rename_map["peso"] = "peso_kg"
    if "cintura" in base.columns: rename_map["cintura"] = "cintura_cm"
    base = base.rename(columns=rename_map)

    if "altura_cm" not in base.columns and "talla_m" in base.columns:
        base["talla_m"] = base["talla_m"].map(_to_num)
        base["altura_cm"] = base["talla_m"] * 100.0

    required = {"encuesta_id","fotos_frontal","fotos_lateral","peso_kg","altura_cm","cintura_cm"}
    miss = required - set(base.columns)
    if miss:
        raise ValueError(f"Faltan columnas: {miss}")

    for col in ["peso_kg","altura_cm","cintura_cm"]:
        base[col] = base[col].map(_to_num)

    tab_cols = [c.strip() for c in args.tab_cols.split(",")] if args.tab_cols else []
    tab_cols = [c for c in tab_cols if c]

    forbidden = {"peso","peso_kg","altura_cm","talla_m","cintura","cintura_cm","ict","imc","whtr"}
    leak = sorted(set(tab_cols) & forbidden)
    if leak:
        print_rank0(f"Advertencia: eliminando columnas objetivo: {leak}")
        tab_cols = [c for c in tab_cols if c not in forbidden]

    present_tab = []
    for c in tab_cols:
        if c in base.columns:
            if base[c].dtype == object:
                tmp = base[c].map(_to_num)
                na_ratio = tmp.isna().mean()
                if na_ratio <= 0.3:
                    base[c] = tmp
                else:
                    try:
                        base[c] = base[c].astype(float)
                    except Exception:
                        pass
            present_tab.append(c)
        else:
            print_rank0(f"Columna '{c}' no está en CSV. Se ignora.")

    tab_cols = present_tab

    base = base.dropna(subset=["encuesta_id","fotos_frontal","fotos_lateral","peso_kg","altura_cm","cintura_cm"]).copy()
    base = base[(base["peso_kg"]>0) & (base["altura_cm"]>0) & (base["cintura_cm"]>0)].copy()

    df_all = expand_all_photos(base)
    if df_all.empty:
        raise RuntimeError("Sin pares válidos")

    for c in tab_cols:
        df_all[c] = base.set_index("encuesta_id").loc[df_all["encuesta_id"].values, c].values

    groups = df_all["encuesta_id"].astype(str).values
    splitter = GroupShuffleSplit(n_splits=1, train_size=0.8, random_state=args.seed)
    train_idx, val_idx = next(splitter.split(df_all, groups=groups))
    df_train = df_all.iloc[train_idx].reset_index(drop=True)
    df_val   = df_all.iloc[val_idx].reset_index(drop=True)

    tab_medians = {}
    if len(tab_cols):
        for c in tab_cols:
            df_train[c] = pd.to_numeric(df_train[c], errors="coerce")
            df_val[c]   = pd.to_numeric(df_val[c], errors="coerce")
            med = df_train[c].median()
            tab_medians[c] = float(med) if np.isfinite(med) else 0.0
            df_train[c] = df_train[c].fillna(tab_medians[c])
            df_val[c]   = df_val[c].fillna(tab_medians[c])

    tmp_ds = DuoViewsImgTabDataset(df_train, tab_cols, train=True, img_size=args.img_size)
    y_mean, y_std = tmp_ds.y_mean, tmp_ds.y_std
    tab_mean, tab_std = tmp_ds.tab_mean, tmp_ds.tab_std
    train_ds = tmp_ds
    val_ds = DuoViewsImgTabDataset(df_val, tab_cols, train=False, img_size=args.img_size,
                                   y_mean=y_mean, y_std=y_std, tab_mean=tab_mean, tab_std=tab_std)

    # IMPORTANTE: DistributedSampler para DDP
    train_sampler = DistributedSampler(train_ds, num_replicas=world_size, rank=rank, shuffle=True)

    train_loader = DataLoader(train_ds, batch_size=args.batch_size, sampler=train_sampler,
                              num_workers=args.workers, pin_memory=True, drop_last=True)
    val_loader = DataLoader(val_ds, batch_size=args.batch_size, shuffle=False,
                            num_workers=args.workers, pin_memory=True)

    print_rank0(f"[Rank {rank}] Train: {len(train_ds)}, Val: {len(val_ds)}")

    # Modelo con DDP
    model = TwoBranchImgTabRegressor(backbone=args.backbone, dropout=args.dropout,
                                     out_dim=3, tab_dim=len(tab_cols), tab_hidden=args.tab_hidden).to(device)
    model = DDP(model, device_ids=[local_rank])

    criterion = nn.SmoothL1Loss(beta=0.5)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-4)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)

    best = 1e9
    for epoch in range(1, args.epochs+1):
        train_sampler.set_epoch(epoch)  # CRÍTICO para DDP

        model.train()
        if rank == 0:
            pbar = tqdm(train_loader, desc=f"Epoch {epoch}/{args.epochs}")
        else:
            pbar = train_loader

        for batch in pbar:
            f = batch["front"].to(device, non_blocking=True)
            s = batch["side"].to(device, non_blocking=True)
            t = batch["tab"].to(device, non_blocking=True) if len(tab_cols) else None
            y = batch["y"].to(device, non_blocking=True)

            optimizer.zero_grad(set_to_none=True)
            yhat = model(f, s, t)

            loss_all = criterion(yhat, y)
            if args.w_ict != 1.0 or args.w_peso != 1.0 or args.w_altura != 1.0:
                l_peso   = criterion(yhat[:,0:1], y[:,0:1])
                l_altura = criterion(yhat[:,1:2], y[:,1:2])
                l_ict    = criterion(yhat[:,2:3], y[:,2:3])
                loss_all = args.w_peso*l_peso + args.w_altura*l_altura + args.w_ict*l_ict

            loss_all.backward()
            nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()

            if rank == 0:
                pbar.set_postfix(loss=f"{loss_all.item():.4f}")

        # Validación (solo rank 0)
        if rank == 0:
            model.eval()
            val_loss = 0.0
            preds, tgts = [], []
            with torch.no_grad():
                for batch in val_loader:
                    f = batch["front"].to(device, non_blocking=True)
                    s = batch["side"].to(device, non_blocking=True)
                    t = batch["tab"].to(device, non_blocking=True) if len(tab_cols) else None
                    y = batch["y"].to(device, non_blocking=True)
                    yhat = model(f, s, t)

                    if args.w_ict != 1.0 or args.w_peso != 1.0 or args.w_altura != 1.0:
                        l_peso   = criterion(yhat[:,0:1], y[:,0:1])
                        l_altura = criterion(yhat[:,1:2], y[:,1:2])
                        l_ict    = criterion(yhat[:,2:3], y[:,2:3])
                        loss = args.w_peso*l_peso + args.w_altura*l_altura + args.w_ict*l_ict
                    else:
                        loss = criterion(yhat, y)

                    val_loss += loss.item() * y.size(0)
                    preds.append(yhat.cpu().numpy()); tgts.append(y.cpu().numpy())

            val_loss /= max(1, len(val_ds))

            preds = np.concatenate(preds, axis=0) * val_ds.y_std + val_ds.y_mean
            tgts  = np.concatenate(tgts,  axis=0) * val_ds.y_std + val_ds.y_mean

            mae_w   = mae(preds[:,0], tgts[:,0])
            rmse_w  = rmse(preds[:,0], tgts[:,0])
            mae_h   = mae(preds[:,1], tgts[:,1])
            rmse_h  = rmse(preds[:,1], tgts[:,1])
            mae_ict = mae(preds[:,2], tgts[:,2])
            rmse_ict= rmse(preds[:,2], tgts[:,2])

            print(f"\nEpoch {epoch} Val: loss={val_loss:.4f} | "
                  f"MAE peso={mae_w:.2f}kg RMSE={rmse_w:.2f} | "
                  f"MAE altura={mae_h:.2f}cm RMSE={rmse_h:.2f} | "
                  f"MAE ICT={mae_ict:.4f} RMSE={rmse_ict:.4f}")

            if val_loss < best:
                best = val_loss
                ckpt = {
                    "state_dict": model.module.state_dict(),  # .module para DDP
                    "backbone": args.backbone,
                    "img_size": args.img_size,
                    "y_mean": val_ds.y_mean,
                    "y_std":  val_ds.y_std,
                    "tab_cols": tab_cols,
                    "tab_mean": val_ds.tab_mean if val_ds.tab_mean is None else val_ds.tab_mean.squeeze().tolist(),
                    "tab_std":  val_ds.tab_std  if val_ds.tab_std  is None else val_ds.tab_std.squeeze().tolist(),
                    "tab_medians": tab_medians,
                    "targets": ["peso_kg","altura_cm","ict"],
                }
                path = os.path.join(args.out_dir, "best_ddp_model.pt")
                torch.save(ckpt, path)
                print(f"✅ Guardado: {path}")

        scheduler.step()
        dist.barrier()  # Sincronizar todos los ranks

    # Guardar métricas finales (solo rank 0)
    if rank == 0:
        metrics = {
            "val_loss_norm": float(val_loss),
            "MAE": {"peso_kg": mae_w, "altura_cm": mae_h, "ict": mae_ict},
            "RMSE": {"peso_kg": rmse_w, "altura_cm": rmse_h, "ict": rmse_ict},
            "train_size": int(len(train_ds)), "val_size": int(len(val_ds)),
        }
        with open(os.path.join(args.out_dir, "metrics_ddp.json"), "w") as f:
            json.dump(metrics, f, indent=2)
        print("🎉 Entrenamiento DDP completado!")

    dist.destroy_process_group()

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv_path", default="/mnt/data/dataset_both_views.csv")
    ap.add_argument("--out_dir",  default="./runs_ddp")
    ap.add_argument("--img_size", type=int, default=224)
    ap.add_argument("--batch_size", type=int, default=32)
    ap.add_argument("--epochs", type=int, default=20)
    ap.add_argument("--lr", type=float, default=2e-4)
    ap.add_argument("--backbone", choices=["efficientnet_b0","efficientnet_b2","efficientnet_b3","resnet18"], default="resnet18")
    ap.add_argument("--dropout", type=float, default=0.3)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--tab_cols", type=str, default="")
    ap.add_argument("--tab_hidden", type=int, default=64)
    ap.add_argument("--w_peso", type=float, default=1.0)
    ap.add_argument("--w_altura", type=float, default=1.0)
    ap.add_argument("--w_ict", type=float, default=1.5)
    args = ap.parse_args()
    main(args)
