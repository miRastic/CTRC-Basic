# conda create -n Exologer2 python=3.8.14
# conda run -n Exologer2 python ./ExoLoger-v1.py
# conda run -n Exologer2 python /mnt/DataDrive/MDD/T0039--PROJ005--Exo-miR-DB-Py/ExoLoger-Server/ExoLoger-v1.py
# requirements: conda install UtilityLib mechanicalsoup tensorflow sklearn

from UtilityLib import ProjectManager
from SeqUtility import RNASeqEncoder
# import tensorflow as TF
# import keras.models as K_MODELS

Enc_rna = RNASeqEncoder()
EXOLOG2 = ProjectManager(
  path_bases=["/mnt/SuperData/ExoLoger-Server",
              r"U:\MDD-CTRL\T0039--PROJ005--Exo-miR-DB-Py\ExoLoger-Server"],
  subversion="server24")

EXOLOG2.require_many([
    ("pandas", "PD"), ("json", "JSON"), ('math', "MATH"), ('mechanicalsoup', 'MS'),
    ("tensorflow", "TF"), ("keras.models", "K_MODELS")
  ])

_model_id = "V5-D16G97e18"

# @TF.function
def predict_is_exo(_seq, processed=False):
  global _model
  if not processed:
    _seq = process_seq(_seq)
  _result = _model.predict(_seq)
  return _result[0]

# @TF.function
def process_seq(_seq):
  global _required_shape
  _seq = Enc_rna.encode_seq(_seq[:_required_shape], fill_len=_required_shape)
  _seq = EXOLOG2.TF.convert_to_tensor([_seq])
  return _seq

_browser = EXOLOG2.MS.StatefulBrowser(user_agent='Mozilla/5.0 (platform; rv:geckoversion) Gecko/geckotrail AUTOJOB/VIVIS')

base_url = "https://www.mirna.in/api/server/exologer"
# base_url = "http://ctrl.vishalkumarsahu.in/api/server/exologer"

jobs_url = f"{base_url}/pygetjb23273908/html"
update_url = f"{base_url}/supdjbs23273908/%s"
update_url_status = f"{update_url}/%s"

_table_raw = _browser.get(jobs_url)

_table_list = EXOLOG2.PD.read_html(_table_raw.text)
_jobs = _table_list[0] if len(_table_list) > 0 else None

print(f"There are {_jobs.shape[0]} new job(s).")

for _, _job in EXOLOG2.ProgressBar(_jobs.iterrows(), total=_jobs.shape[0]):
  _res = _browser.get(update_url_status % (_job.job_id, "2"))

EXOLOG2.config.jobs_done = EXOLOG2.config.jobs_done if isinstance(EXOLOG2.config.jobs_done, (int)) else 0
EXOLOG2.update_config()

_model = EXOLOG2.K_MODELS.load_model(EXOLOG2.get_path("Models/V5-D16G97e18.h5"))
_required_shape = _model.layers[0].input_shape[1]

__loop = _jobs.copy()
for _, _job in EXOLOG2.ProgressBar(__loop.iterrows(), total=__loop.shape[0]):
  _ie, _ne = predict_is_exo(_job.seq, False)

  _confidence = round(sum([_ie, _ne])*100, 2)
  _confusion = round(_ne*100, 2)
  _exosomal = round(_ie*100, 3)
  _is_exo = 1 if (_ie > _ne) else 0
  _res = _browser.post(update_url % 'x-s3cr3t',
      headers=_browser.session.cookies,
      data={
        "job_id": _job.job_id,
        "status": "4",
        "prediction": EXOLOG2.JSON.dumps({
          "model": _model_id,
          "exo_percentage": _exosomal,
          "is_exosomal": _is_exo,
          "confidence": _confidence,
          # "confusion": _confusion,
        }),
      })
  # print(_job.seq, _ie, _ne)
  # print(_res)
  EXOLOG2.config.jobs_done += 1
  EXOLOG2.update_config()
  EXOLOG2.time_sleep(0.1)
