"""
_bigbed_reader.py — Pure-Python BigBed reader (no external tools required)
===========================================================================
Implements a minimal BigBed parser sufficient to extract all records from
ENCODE cCRE V3 and other BigBed files.

BigBed format reference:
  https://genome.ucsc.edu/goldenPath/help/bigBed.html
  https://github.com/ucscGenomeBrowser/kent/blob/master/src/lib/bbiFile.c

Format summary
--------------
File layout (offsets from file start):
  0        : Header (64 bytes)
  chromTreeOffset : B+ tree of chromosome names (chrom → id, size)
  dataOffset      : Compressed data blocks (zlib)
  indexOffset     : R tree index of data blocks

Header (little-endian):
  uint32  magic           = 0x8789F2EB
  uint16  version         = 4
  uint16  zoomLevels
  uint64  chromTreeOffset
  uint64  unzoomedDataOffset
  uint64  unzoomedIndexOffset
  uint16  fieldCount
  uint16  definedFieldCount
  uint64  autoSqlOffset
  uint64  totalSummaryOffset
  uint32  uncompressBufSize   (0 = uncompressed)
  uint64  extensionOffset

B+ tree node:
  uint8   isLeaf
  uint8   reserved
  uint16  count
  [Leaf  items]: key[keySize] + uint32 chromId + uint32 chromSize
  [Inner items]: key[keySize] + uint64 childOffset

R tree node:
  uint8   isLeaf
  uint8   reserved
  uint16  count
  [Leaf  items]: startChromIx(4) startBase(4) endChromIx(4) endBase(4)
                 dataOffset(8) dataSize(8)  = 32 bytes
  [Inner items]: startChromIx(4) startBase(4) endChromIx(4) endBase(4)
                 childOffset(8)              = 24 bytes

Data records (after zlib decompression):
  uint32  chromId
  uint32  chromStart
  uint32  chromEnd
  char[]  rest   (null-terminated; extra fields tab-separated)
"""
import struct
import zlib
from typing import Dict, Generator, List, Optional, Tuple

# ── Constants ─────────────────────────────────────────────────────────────────
BIGBED_MAGIC    = 0x8789F2EB
CHROM_TREE_MAGIC = 0x78CA8C91
R_TREE_MAGIC    = 0x2468ACE0


# ── Header ────────────────────────────────────────────────────────────────────

def _read_header(f) -> dict:
    raw = f.read(64)
    if len(raw) < 64:
        raise ValueError("File too short to be a BigBed")
    magic, version, zoom_levels = struct.unpack_from('<IHH', raw, 0)
    if magic != BIGBED_MAGIC:
        raise ValueError(f"Not a BigBed file (magic={hex(magic)})")
    chrom_tree_off, data_off, index_off = struct.unpack_from('<QQQ', raw, 8)
    field_count, defined_field_count    = struct.unpack_from('<HH', raw, 32)
    autosql_off, total_summary_off      = struct.unpack_from('<QQ', raw, 36)
    uncompress_buf_size, extension_off  = struct.unpack_from('<IQ', raw, 52)
    return {
        "version":            version,
        "zoom_levels":        zoom_levels,
        "chrom_tree_off":     chrom_tree_off,
        "data_off":           data_off,
        "index_off":          index_off,
        "field_count":        field_count,
        "defined_field_count": defined_field_count,
        "autosql_off":        autosql_off,
        "uncompress_buf_size": uncompress_buf_size,
        "extension_off":      extension_off,
    }


# ── B+ tree (chromosome name lookup) ─────────────────────────────────────────

def _read_chrom_tree(f, tree_offset: int) -> Dict[int, str]:
    """Return {chromId: chrom_name} mapping from the B+ tree."""
    f.seek(tree_offset)
    magic, block_size, key_size, val_size, item_count, reserved = struct.unpack('<IIIIQQ', f.read(32))
    if magic != CHROM_TREE_MAGIC:
        raise ValueError(f"Bad chromTree magic: {hex(magic)}")

    chroms: Dict[int, str] = {}
    _read_chrom_node(f, key_size, val_size, chroms)
    return chroms


def _read_chrom_node(f, key_size: int, val_size: int, chroms: Dict[int, str]) -> None:
    is_leaf, _reserved, count = struct.unpack('<BBH', f.read(4))
    if is_leaf:
        for _ in range(count):
            key_bytes = f.read(key_size)
            chrom_id, chrom_size = struct.unpack('<II', f.read(val_size))
            chrom_name = key_bytes.rstrip(b'\x00').decode('ascii', errors='replace')
            chroms[chrom_id] = chrom_name
    else:
        child_offsets = []
        for _ in range(count):
            _key_bytes = f.read(key_size)
            child_off = struct.unpack('<Q', f.read(8))[0]
            child_offsets.append(child_off)
        for off in child_offsets:
            f.seek(off)
            _read_chrom_node(f, key_size, val_size, chroms)


# ── R tree (data block index) ─────────────────────────────────────────────────

def _collect_leaf_blocks(f, node_offset: int) -> List[Tuple[int, int]]:
    """
    Recursively traverse the R tree, returning list of (dataOffset, dataSize)
    for every leaf node.
    """
    f.seek(node_offset)
    is_leaf, _reserved, count = struct.unpack('<BBH', f.read(4))
    blocks: List[Tuple[int, int]] = []

    if is_leaf:
        for _ in range(count):
            # startChromIx(4) startBase(4) endChromIx(4) endBase(4) dataOffset(8) dataSize(8)
            _start_chrom, _start_base, _end_chrom, _end_base, data_off, data_size = \
                struct.unpack('<IIIIQQ', f.read(32))
            blocks.append((data_off, data_size))
    else:
        child_offsets: List[int] = []
        for _ in range(count):
            # startChromIx(4) startBase(4) endChromIx(4) endBase(4) childOffset(8)
            _sc, _sb, _ec, _eb, child_off = struct.unpack('<IIIIq', f.read(24))
            child_offsets.append(child_off)
        # Save position before recursing
        for child_off in child_offsets:
            blocks.extend(_collect_leaf_blocks(f, child_off))

    return blocks


def _read_r_tree_header(f, tree_offset: int) -> Tuple[int, List[Tuple[int, int]]]:
    """
    Read R tree header and collect all leaf data blocks.
    Returns (items_per_slot, [(dataOffset, dataSize), ...])
    """
    f.seek(tree_offset)
    magic, block_size, item_count = struct.unpack('<IIQ', f.read(16))
    if magic != R_TREE_MAGIC:
        raise ValueError(f"Bad rTree magic: {hex(magic)}")
    # Skip: startChromIx(4) startBase(4) endChromIx(4) endBase(4) endFileOffset(8)
    f.read(24)
    items_per_slot, _reserved = struct.unpack('<II', f.read(8))

    root_offset = f.tell()  # root node starts right after the header
    blocks = _collect_leaf_blocks(f, root_offset)
    return items_per_slot, blocks


# ── Data block parser ─────────────────────────────────────────────────────────

def _parse_block(
    block_data: bytes,
    chroms: Dict[int, str],
) -> Generator[Tuple[str, int, int, List[str]], None, None]:
    """
    Parse a (decompressed) data block.

    Each record:
      uint32  chromId
      uint32  chromStart
      uint32  chromEnd
      char[]  rest  (null-terminated string, fields tab-separated)

    Yields: (chrom_name, start, end, extra_fields_list)
    """
    pos = 0
    n = len(block_data)
    while pos < n:
        if pos + 12 > n:
            break
        chrom_id, chrom_start, chrom_end = struct.unpack_from('<III', block_data, pos)
        pos += 12

        # Read null-terminated string for extra fields
        end = block_data.find(b'\x00', pos)
        if end == -1:
            rest_bytes = block_data[pos:]
            pos = n
        else:
            rest_bytes = block_data[pos:end]
            pos = end + 1

        chrom_name = chroms.get(chrom_id, f"chr{chrom_id}")
        extra = rest_bytes.decode('utf-8', errors='replace').split('\t')
        yield chrom_name, chrom_start, chrom_end, extra


# ── Public API ────────────────────────────────────────────────────────────────

def iter_records(
    path: str,
) -> Generator[Tuple[str, int, int, List[str]], None, None]:
    """
    Iterate over all records in a BigBed file.

    Yields: (chrom_name, chromStart, chromEnd, extra_fields)
      chrom_name  : e.g. "chr1"
      chromStart  : 0-based
      chromEnd    : 0-based exclusive
      extra_fields: list of strings (fields after chrom/start/end)

    Records are yielded in no guaranteed order (block order from R tree).
    """
    with open(path, 'rb') as f:
        hdr = _read_header(f)
        chroms = _read_chrom_tree(f, hdr["chrom_tree_off"])
        _items_per_slot, blocks = _read_r_tree_header(f, hdr["index_off"])
        do_decompress = hdr["uncompress_buf_size"] > 0

        for data_off, data_size in blocks:
            f.seek(data_off)
            raw = f.read(data_size)
            if do_decompress:
                try:
                    raw = zlib.decompress(raw)
                except zlib.error as e:
                    continue  # skip corrupt block
            yield from _parse_block(raw, chroms)


def field_count(path: str) -> Tuple[int, int]:
    """Return (fieldCount, definedFieldCount) from the BigBed header."""
    with open(path, 'rb') as f:
        hdr = _read_header(f)
    return hdr["field_count"], hdr["defined_field_count"]


def autosql(path: str) -> Optional[str]:
    """Return the autoSql schema string, or None if not present."""
    with open(path, 'rb') as f:
        hdr = _read_header(f)
        if hdr["autosql_off"] == 0:
            return None
        f.seek(hdr["autosql_off"])
        raw = f.read(4096)
        nul = raw.find(b'\x00')
        if nul != -1:
            raw = raw[:nul]
        return raw.decode('utf-8', errors='replace')
