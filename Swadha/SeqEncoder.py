import numpy
from sklearn.preprocessing import OneHotEncoder
from UtilityLib import EU
class SeqEncoder():
  def __init__(self, *args, **kwargs):
    self.__defaults = {
        "debug": False,
        "model": None,
        "NP": numpy,
        "residues": "AUTGCX", # AUTGCX
        "residue_unknown": "X",
      }
    EU.update_attributes(self, kwargs, self.__defaults)

  def __call__(self, *args, **kwargs):
    return self.encode_seq(*args, **kwargs)

  def process_seq(self, *args, **kwargs):
    _seq = args[0] if len(args) > 0 else kwargs.get("seq")
    _fill_len = args[1] if len(args) > 1 else kwargs.get("fill_len")
    _fill_align = args[2] if len(args) > 2 else kwargs.get("fill_align", "left")

    if _fill_len:
      if _fill_align == "left":
        _seq = _seq.ljust(_fill_len, self.residue_unknown)
      elif _fill_align == "center":
        _seq = _seq.center(_fill_len, self.residue_unknown)
      else:
        _seq = _seq.rjust(_fill_len, self.residue_unknown)

    # Validate sequence
    _seq = [_s if _s in self.residues else self.residue_unknown for _s in _seq]
    _seq = self.NP.array(_seq).reshape(-1, 1)
    return _seq

  def init_model(self, *args, **kwargs):
    _seq_restricted = self.process_seq(self.residues, **kwargs)
    self.model = OneHotEncoder()
    self.model.fit_transform(_seq_restricted)

  def encode_seq(self, *args, **kwargs):
    _seq = args[0] if len(args) > 0 else kwargs.get("seq")
    _seq = self.process_seq(_seq, **kwargs)
    # https://stackoverflow.com/a/68689457
    # return self.NP.array([self.NP.array(_v) for _v in self.model.transform(_seq).toarray()])

    _seq = self.model.transform(_seq).toarray()
    # _enc = []
    # for _s in _seq:
    #   _enc.extend(_s)
    return _seq

# RNA Secondary Structure Seq Encoder
class RNASSEncoder(SeqEncoder):
  def __init__(self, *args, **kwargs):
    super(RNASSEncoder, self).__init__(**kwargs)
    kwargs["residues"] = ".()}{[]" #
    kwargs["residue_unknown"] = "_" #
    self.__defaults = {}
    EU.update_attributes(self, kwargs, self.__defaults)
    self.init_model()

# Nucleic Acid Seq Encoder
class NASeqEncoder(SeqEncoder):
  def __init__(self, *args, **kwargs):
    super(NASeqEncoder, self).__init__(**kwargs)
    kwargs["residues"] = "AUTGCX"
    self.__defaults = {}
    EU.update_attributes(self, kwargs, self.__defaults)
    self.init_model()

# RNA Seq Encoder
class RNASeqEncoder(NASeqEncoder):
  def __init__(self, *args, **kwargs):
    super(RNASeqEncoder, self).__init__(**kwargs)
    kwargs["residues"] = "AUGCX"
    self.__defaults = {}
    EU.update_attributes(self, kwargs, self.__defaults)
    self.init_model()

# DNA Seq Encoder
class DNASeqEncoder(NASeqEncoder):
  def __init__(self, *args, **kwargs):
    super(DNASeqEncoder, self).__init__(**kwargs)
    kwargs["residues"] = "ATGCX"
    self.__defaults = {}
    EU.update_attributes(self, kwargs, self.__defaults)
    self.init_model()

# Protein/peptide/AA Seq Encoder
class AASeqEncoder(SeqEncoder):
  def __init__(self, *args, **kwargs):
    super(AASeqEncoder, self).__init__(**kwargs)
    # AA Except BJOUWXZ (X represents other)
    kwargs["residues"] = "ABCDEFGHIKLMNPQRSTUVWXYZ"
    self.__defaults = {}
    EU.update_attributes(self, kwargs, self.__defaults)
    self.init_model()
