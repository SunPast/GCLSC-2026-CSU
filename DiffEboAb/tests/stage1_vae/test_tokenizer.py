"""Unit tests for AntibodyTokenizer."""
from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.data.tokenizer import AntibodyTokenizer, build_tokenizer


def test_vocab_size_and_special_ids():
    tok = build_tokenizer(max_len=140)
    assert tok.vocab_size == 26
    assert tok.pad_id == 0
    assert tok.n_special == 6
    assert set(tok.AMINO_ACIDS) <= set(tok._stoi)


def test_encode_decode_roundtrip_aa():
    tok = build_tokenizer()
    seq = "EVQLVESGGGLVQPGGSLRLSCAAS"
    ids = tok.encode(seq, chain="H")
    assert ids[0] == tok.h_id
    assert ids[-1] == tok.eos_id
    assert tok.decode(ids) == seq


def test_encode_light_chain_tag():
    tok = build_tokenizer()
    ids = tok.encode("DIQMTQSPSSLSASVGDR", chain="L")
    assert ids[0] == tok.l_id


def test_decode_stops_at_eos():
    tok = build_tokenizer()
    ids = [tok.h_id, tok.token_to_id("E"), tok.token_to_id("V"), tok.eos_id, tok.token_to_id("Q")]
    assert tok.decode(ids, stop_at_eos=True) == "EV"


def test_truncate_keeps_eos():
    tok = AntibodyTokenizer(max_len=8)
    ids = tok.encode("ACDEFGHIKLMN", chain="H", max_len=8)
    assert len(ids) == 8
    assert ids[-1] == tok.eos_id


def test_is_valid_aa():
    tok = build_tokenizer()
    assert tok.is_valid_aa("ACDE")
    assert not tok.is_valid_aa("ACDEX")  # X not in alphabet as AA letter for is_valid
