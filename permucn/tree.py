"""Phylogeny parsing and branch-index utilities."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


_STATE_SUFFIX = re.compile(r"^(.*)_([0-9]+)$")
_TIP_SUFFIX = re.compile(r"^(.*)<[0-9]+>$")


@dataclass(frozen=True)
class CanonicalTree:
    root: int
    labels: List[str]
    branch_keys_by_node: List[str]
    parent_by_node: List[Optional[int]]
    children_by_node: List[List[int]]
    branch_length_by_node: List[float]
    tip_species_by_node: List[Optional[str]]
    non_root_nodes: List[int]
    branch_key_by_index: List[str]
    node_by_branch_index: List[int]
    branch_index_by_key: Dict[str, int]
    desc_mask_by_branch_index: List[int]
    anc_mask_by_branch_index: List[int]
    clade_size_by_branch_index: List[int]
    clade_bin_by_branch_index: List[int]
    all_mask: int


@dataclass
class _Node:
    label: str
    length: float
    children: List[int]


class _NewickParser:
    def __init__(self, text: str) -> None:
        self.text = text
        self.i = 0
        self.n = len(text)
        self.nodes: List[_Node] = []

    def parse(self) -> Tuple[int, List[_Node]]:
        root = self._parse_subtree()
        self._skip_ws()
        if self.i < self.n and self.text[self.i] == ";":
            self.i += 1
        self._skip_ws()
        if self.i != self.n:
            raise ValueError(f"Unexpected trailing text in Newick at pos {self.i}")
        return root, self.nodes

    def _peek(self) -> str:
        if self.i >= self.n:
            return ""
        return self.text[self.i]

    def _consume(self, ch: str) -> None:
        if self._peek() != ch:
            raise ValueError(f"Expected '{ch}' at pos {self.i}, found '{self._peek()}'")
        self.i += 1

    def _skip_ws(self) -> None:
        while self.i < self.n and self.text[self.i].isspace():
            self.i += 1

    def _parse_label(self) -> str:
        self._skip_ws()
        start = self.i
        while self.i < self.n and self.text[self.i] not in ",():;":
            self.i += 1
        return self.text[start:self.i].strip()

    def _parse_length(self) -> float:
        self._skip_ws()
        if self._peek() != ":":
            return float("nan")
        self.i += 1
        self._skip_ws()
        start = self.i
        while self.i < self.n and self.text[self.i] not in ",();":
            self.i += 1
        raw = self.text[start:self.i].strip()
        if not raw:
            return float("nan")
        return float(raw)

    def _parse_subtree(self) -> int:
        self._skip_ws()
        if self._peek() == "(":
            self._consume("(")
            children: List[int] = []
            while True:
                child = self._parse_subtree()
                children.append(child)
                self._skip_ws()
                if self._peek() == ",":
                    self.i += 1
                    continue
                break
            self._skip_ws()
            self._consume(")")
            label = self._parse_label()
            length = self._parse_length()
            node_id = len(self.nodes)
            self.nodes.append(_Node(label=label, length=length, children=children))
            return node_id

        label = self._parse_label()
        if not label:
            raise ValueError(f"Empty tip label at pos {self.i}")
        length = self._parse_length()
        node_id = len(self.nodes)
        self.nodes.append(_Node(label=label, length=length, children=[]))
        return node_id


def branch_key_from_label(label: str) -> str:
    """Strip CAFE state suffix from a node label.

    Example: "Acanthochromis_polyacanthus<66>_1" -> "Acanthochromis_polyacanthus<66>"
    """
    label = label.strip()
    m = _STATE_SUFFIX.match(label)
    if m and m.group(1):
        return m.group(1)
    return label


def tip_species_from_branch_key(branch_key: str) -> Optional[str]:
    """Extract species name from a tip-like branch key.

    Example: "Acanthochromis_polyacanthus<66>" -> "Acanthochromis_polyacanthus".
    Returns None for internal-like keys such as "<12>".
    """
    m = _TIP_SUFFIX.match(branch_key)
    if m:
        return m.group(1)
    if branch_key.startswith("<") and branch_key.endswith(">"):
        return None
    return branch_key or None


def log2_clade_bin(size: int) -> int:
    if size <= 0:
        raise ValueError(f"Invalid clade size: {size}")
    if size == 1:
        return 0
    if size == 2:
        return 1
    if size <= 4:
        return 2
    if size <= 8:
        return 3
    if size <= 16:
        return 4
    if size <= 32:
        return 5
    if size <= 64:
        return 6
    return 7


def read_first_tree_newick(path: str | Path) -> str:
    """Extract first TREE entry from a NEXUS file."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Tree file not found: {path}")

    started = False
    chunks: List[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not started:
                if "TREE" in line and "=" in line:
                    started = True
                    after = line.split("=", 1)[1].strip()
                    chunks.append(after)
                    if ";" in after:
                        break
            else:
                part = line.strip()
                chunks.append(part)
                if ";" in part:
                    break

    if not chunks:
        raise ValueError(f"No TREE entry found in: {path}")

    text = " ".join(chunks)
    semi = text.find(";")
    if semi == -1:
        raise ValueError(f"Malformed TREE entry (missing ';') in: {path}")
    return text[: semi + 1]


def load_canonical_tree(path: str | Path) -> CanonicalTree:
    """Parse the first ASR tree and build branch-index caches."""
    newick = read_first_tree_newick(path)
    parser = _NewickParser(newick)
    root, nodes = parser.parse()

    n = len(nodes)
    parent_by_node: List[Optional[int]] = [None] * n
    children_by_node: List[List[int]] = [[] for _ in range(n)]
    labels: List[str] = [""] * n
    branch_keys_by_node: List[str] = [""] * n
    branch_length_by_node: List[float] = [float("nan")] * n
    tip_species_by_node: List[Optional[str]] = [None] * n

    for node_id, node in enumerate(nodes):
        labels[node_id] = node.label
        branch_keys_by_node[node_id] = branch_key_from_label(node.label)
        children_by_node[node_id] = list(node.children)
        branch_length_by_node[node_id] = node.length
        for child in node.children:
            parent_by_node[child] = node_id

    for node_id, children in enumerate(children_by_node):
        if not children:
            tip_species_by_node[node_id] = tip_species_from_branch_key(branch_keys_by_node[node_id])

    non_root_nodes = [i for i in range(n) if i != root]
    node_by_branch_index: List[int] = []
    branch_key_by_index: List[str] = []
    branch_index_by_key: Dict[str, int] = {}

    for node_id in non_root_nodes:
        key = branch_keys_by_node[node_id]
        if not key:
            raise ValueError(f"Encountered empty branch key on non-root node id {node_id}")
        if key in branch_index_by_key:
            raise ValueError(f"Duplicate branch key in tree: {key}")
        idx = len(node_by_branch_index)
        node_by_branch_index.append(node_id)
        branch_key_by_index.append(key)
        branch_index_by_key[key] = idx

    m = len(node_by_branch_index)

    node_to_branch_idx: Dict[int, int] = {node_id: idx for idx, node_id in enumerate(node_by_branch_index)}

    tip_count_by_node: List[int] = [0] * n
    desc_mask_by_node: List[int] = [0] * n

    postorder = _postorder(root, children_by_node)
    for node_id in postorder:
        children = children_by_node[node_id]
        if not children:
            tip_count_by_node[node_id] = 1
        else:
            tip_count_by_node[node_id] = sum(tip_count_by_node[ch] for ch in children)

        mask = 0
        if node_id != root:
            mask |= 1 << node_to_branch_idx[node_id]
        for ch in children:
            mask |= desc_mask_by_node[ch]
        desc_mask_by_node[node_id] = mask

    desc_mask_by_branch_index = [0] * m
    clade_size_by_branch_index = [0] * m
    clade_bin_by_branch_index = [0] * m
    for idx, node_id in enumerate(node_by_branch_index):
        desc_mask_by_branch_index[idx] = desc_mask_by_node[node_id]
        clade_size = tip_count_by_node[node_id]
        clade_size_by_branch_index[idx] = clade_size
        clade_bin_by_branch_index[idx] = log2_clade_bin(clade_size)

    anc_mask_by_branch_index = [0] * m
    for idx, node_id in enumerate(node_by_branch_index):
        p = parent_by_node[node_id]
        anc = 0
        while p is not None:
            if p != root:
                anc_idx = node_to_branch_idx[p]
                anc |= 1 << anc_idx
            p = parent_by_node[p]
        anc_mask_by_branch_index[idx] = anc

    all_mask = (1 << m) - 1 if m > 0 else 0

    return CanonicalTree(
        root=root,
        labels=labels,
        branch_keys_by_node=branch_keys_by_node,
        parent_by_node=parent_by_node,
        children_by_node=children_by_node,
        branch_length_by_node=branch_length_by_node,
        tip_species_by_node=tip_species_by_node,
        non_root_nodes=non_root_nodes,
        branch_key_by_index=branch_key_by_index,
        node_by_branch_index=node_by_branch_index,
        branch_index_by_key=branch_index_by_key,
        desc_mask_by_branch_index=desc_mask_by_branch_index,
        anc_mask_by_branch_index=anc_mask_by_branch_index,
        clade_size_by_branch_index=clade_size_by_branch_index,
        clade_bin_by_branch_index=clade_bin_by_branch_index,
        all_mask=all_mask,
    )


def _postorder(root: int, children_by_node: List[List[int]]) -> List[int]:
    out: List[int] = []

    def rec(node_id: int) -> None:
        for ch in children_by_node[node_id]:
            rec(ch)
        out.append(node_id)

    rec(root)
    return out
