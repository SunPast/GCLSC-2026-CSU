"""Amino-acid tokenizer for BCR/antibody sequences.

Design notes (aligned with AntiBARTy empirical choices):
- Vocabulary: 20 canonical amino acids + MASK, plus special tokens.
  X / unknown residues are filtered upstream in the data pipeline.
- Chain tags: <H> / <L> prefix tokens so the model knows chain type.
- Sequence ends with <EOS>, padding uses <PAD>.
- <MASK> is reserved for the VAE's optional denoising/MLM-style aux loss.

Token indices are stable and exposed as class constants so that
checkpoints/configs stay reproducible across runs.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Sequence


class AntibodyTokenizer:
    """Deterministic residue-level tokenizer for antibody Fv sequences."""

    # ---- special tokens (indices are fixed, do not reorder) ----
    PAD = "<pad>"
    MASK = "<mask>"
    UNK = "<unk>"
    EOS = "<eos>"
    H_CHAIN = "<H>"
    L_CHAIN = "<L>"

    # ---- canonical amino acids (single-letter) ----
    AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

    def __init__(self, max_len: int = 140) -> None:
        # special tokens first, then the 20 AA
        self._itos = [self.PAD, self.MASK, self.UNK, self.EOS,
                      self.H_CHAIN, self.L_CHAIN] + list(self.AMINO_ACIDS)
        self._stoi = {s: i for i, s in enumerate(self._itos)}
        self.max_len = max_len

        # convenience aliases
        self.pad_id = self._stoi[self.PAD]
        self.mask_id = self._stoi[self.MASK]
        self.unk_id = self._stoi[self.UNK]
        self.eos_id = self._stoi[self.EOS]
        self.h_id = self._stoi[self.H_CHAIN]
        self.l_id = self._stoi[self.L_CHAIN]

    # ------------------------------------------------------------------
    # vocabulary introspection
    # ------------------------------------------------------------------
    @property
    def vocab_size(self) -> int:
        return len(self._itos)

    @property
    def n_special(self) -> int:
        """Number of non-amino-acid tokens (specials)."""
        return 6

    def id_to_token(self, idx: int) -> str:
        return self._itos[idx]

    def token_to_id(self, token: str) -> int:
        return self._stoi[token]

    def __len__(self) -> int:
        return self.vocab_size

    # ------------------------------------------------------------------
    # encoding
    # ------------------------------------------------------------------
    def encode(self, seq: str, chain: str = "H",
               add_eos: bool = True, max_len: Optional[int] = None) -> List[int]:
        """Tokenize an amino-acid sequence into token ids.

        Args:
            seq: amino-acid sequence (canonical 20-letter alphabet).
            chain: "H" or "L" -> prepends the matching chain tag.
            add_eos: append <EOS> token.
            max_len: hard cap on total length (incl. tag/EOS); truncates.
        """
        max_len = max_len or self.max_len
        chain_id = self.h_id if chain in ("H", "Heavy", "heavy") else self.l_id
        ids = [chain_id]

        for ch in seq:
            if ch in self._stoi:
                ids.append(self._stoi[ch])
            else:
                ids.append(self.unk_id)  # should not happen after X filtering

        if add_eos:
            ids.append(self.eos_id)

        if len(ids) > max_len:
            ids = ids[: max_len]
            if ids[-1] != self.eos_id:
                ids[-1] = self.eos_id
        return ids

    def encode_batch(self, seqs: Sequence[str], chains: Sequence[str],
                     add_eos: bool = True, max_len: Optional[int] = None
                     ) -> List[List[int]]:
        return [self.encode(s, c, add_eos=add_eos, max_len=max_len)
                for s, c in zip(seqs, chains)]

    # ------------------------------------------------------------------
    # decoding
    # ------------------------------------------------------------------
    def decode(self, ids: Iterable[int], skip_special: bool = True,
               stop_at_eos: bool = True) -> str:
        """Convert token ids back to an amino-acid string.

        Args:
            skip_special: drop PAD/MASK/chain-tag tokens from output.
            stop_at_eos: stop decoding at the first <EOS>.
        """
        out: List[str] = []
        for i in ids:
            tok = self._itos[i]
            if stop_at_eos and tok == self.EOS:
                break
            if skip_special and tok in (self.PAD, self.MASK, self.UNK,
                                        self.H_CHAIN, self.L_CHAIN):
                continue
            out.append(tok)
        return "".join(out)

    def decode_batch(self, batch_ids: Sequence[Sequence[int]],
                     skip_special: bool = True, stop_at_eos: bool = True
                     ) -> List[str]:
        return [self.decode(ids, skip_special, stop_at_eos)
                for ids in batch_ids]

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------
    def is_valid_aa(self, seq: str) -> bool:
        """True if every residue is a canonical amino acid."""
        return all(ch in self._stoi and ch not in (self.PAD, self.MASK,
                                                   self.UNK, self.EOS)
                   for ch in seq)

    def to_dict(self) -> Dict[str, object]:
        return {"vocab": self._itos, "max_len": self.max_len}

    @classmethod
    def from_dict(cls, cfg: Dict[str, object]) -> "AntibodyTokenizer":
        raw_max_len = cfg.get("max_len", 140)
        max_len = int(raw_max_len) if raw_max_len is not None else 140
        tok = cls(max_len=max_len)
        vocab = cfg.get("vocab")
        if isinstance(vocab, (list, tuple)) and len(vocab) > 0:
            tok._itos = [str(v) for v in vocab]
            tok._stoi = {s: i for i, s in enumerate(tok._itos)}
            tok.pad_id = tok._stoi[tok.PAD]
            tok.mask_id = tok._stoi[tok.MASK]
            tok.unk_id = tok._stoi[tok.UNK]
            tok.eos_id = tok._stoi[tok.EOS]
            tok.h_id = tok._stoi[tok.H_CHAIN]
            tok.l_id = tok._stoi[tok.L_CHAIN]
        return tok


def build_tokenizer(max_len: int = 140) -> AntibodyTokenizer:
    return AntibodyTokenizer(max_len=max_len)
