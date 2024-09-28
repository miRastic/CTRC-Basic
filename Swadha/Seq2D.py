from UtilityLib import ProjectManager
# from BioScrapper import ChromeManager

class SeqScrap(ProjectManager):
  def __init__(self, *args, **kwargs):
    super(SeqScrap, self).__init__(**kwargs)
    self.__defaults = {
      "path_base": ".",
      "required_libs": [("pandas", "PD"),
        ("urllib.request", "URLREQ"),
        ("numpy", "NP"),
        ("requests", "REQUESTS"),
        ("matplotlib.pyplot", "PLOT"),
        ("seaborn", "SNS")],
      "subversion": "SeqScrap",
      "Session_Clone": None,
    }
    self.update_attributes(self, kwargs, self.__defaults)
    self.require_many(self.required_libs)

  def init_scrap(self, *args, **kwargs):
    self.DRIVER = ChromeManager(delay=1, **kwargs)
    self.CHROME_INSTANCE = self.DRIVER.wd_instance

  def finish_scrap(self):
    self.DRIVER.close()

  def clone_session(self):
    self.Session_Clone = self.REQUESTS.session()

    # Set correct user agent
    _user_agent = self.CHROME_INSTANCE.execute_script("return navigator.userAgent;")
    self.Session_Clone.headers.update({"user-agent": _user_agent})

    for cookie in self.CHROME_INSTANCE.get_cookies():
      self.Session_Clone.cookies.set(cookie['name'], cookie['value'], domain=cookie['domain'])

# RNA 2D Structure
class Seq2D_RNA(SeqScrap):
  def __init__(self, *args, **kwargs):
    super(Seq2D_RNA, self).__init__(**kwargs)
    self.__defaults = {}
    self.update_attributes(self, kwargs, self.__defaults)

class Seq2D_MXFold2(Seq2D_RNA):
  def __init__(self, *args, **kwargs):
    super(Seq2D_MXFold2, self).__init__(**kwargs)
    self.__defaults = {
      "subversion": "mxfold2",
      "project_dir": "mxfold2",
      "url_server": "http://www.dna.bio.keio.ac.jp/mxfold2"
    }
    self.update_attributes(self, kwargs, self.__defaults)
    self.storage_path = self.validate_dir(self.get_path(self.project_dir))

  def submit_download_2d(self, _id, _seq):
    if _id in self.config.mxfold2.rna_collected.keys() or self.exists(f"{self.storage_path}/{_id}.txt"):
        return

    _seq_info = f""">{_id}
{_seq.strip()}
"""
    _url_data = self.DRIVER.get_url(self.url_server)
    _query = self.CHROME_INSTANCE.find_elements(self.DRIVER.wd_by.NAME, "seq")[0]
    _query.click()
    _query.send_keys(_seq_info)

    _submit = self.CHROME_INSTANCE.find_elements(self.DRIVER.wd_by.CSS_SELECTOR, "input[type='submit']")[0]
    _submit.click()

    _pre_el = self.CHROME_INSTANCE.find_elements(self.DRIVER.wd_by.CSS_SELECTOR, 'pre.border')[0]
    _fp = f"{self.storage_path}/{_id}.txt"

    self.write(_fp, _pre_el.text)

    self.config.mxfold2.rna_collected[_id] = {
        self.time_stamp()
    }

    self.update_config()

  def parse_2D(self, _rdata):
    self.re_score = self.re_compile(r'(.*)\s\((.*)\)')

    _ss_score_coll = []

    for _idx, _rna_seq in self.ProgressBar(_rdata):
        _text_rna = self.read_text(f"{self.storage_path}/{_idx}.txt")
        _rna_ss, _rna_score = self.re_score.findall(_text_rna[-1])[0]
        _score_parsed = {
            "rna_acc": _idx,
            "rna_ss": _rna_ss,
            "rna_score": _rna_score,
        }
        _ss_score_coll.append(_score_parsed)

    self.df_score = self.DF(_ss_score_coll)
    self.df_score.to_pickle(self.get_path("mxfold2.df.gz"))
    self.update_config()

  def get_2D(self, _rdata):
    # List of Tuples like ([rna_id, rna_seq], ...)
    self.init_scrap()
    for _id, _seq in self.ProgressBar(_rdata):
      self.submit_download_2d(_id, _seq)
    self.finish_scrap()
    self.parse_2D(_rdata)

class Seq2D_RNAFold(Seq2D_RNA):
  def __init__(self, *args, **kwargs):
    super(Seq2D_RNAFold, self).__init__(**kwargs)
    self.__defaults = {
      "subversion": "rnafold",
      "project_dir": "rnafold",
      "url_server": "http://rna.tbi.univie.ac.at/cgi-bin/RNAWebSuite/RNAfold.cgi"
    }
    self.update_attributes(self, kwargs, self.__defaults)
    self.storage_path = self.validate_dir(self.get_path("rnafold"))

  def submit_2d(self, _id, _seq):
    if _id in self.config.rnafold.rna_submitted.keys() or self.exists(f"{self.storage_path}/{_id}--result.html"):
        return

    _seq_info = f""">{_id}
{_seq.strip()}
"""

    _url_data = self.DRIVER.get_url(self.url_server)

    _query = self.CHROME_INSTANCE.find_elements(self.DRIVER.wd_by.NAME, "SCREEN")[0]
    _query.click()
    _query.send_keys(_seq_info)

    _submit = self.CHROME_INSTANCE.find_elements(self.DRIVER.wd_by.CSS_SELECTOR, "input[type='submit']")[0]
    _submit.click()

    # wait for URL to change with 15 seconds timeout
    self.time_sleep(1)
    # DRIVER.wd_ui.WebDriverWait(self.CHROME_INSTANCE, 15).until(DRIVER.wd_ec.url_changes(_url)) # fix this line

    _anchors = self.CHROME_INSTANCE.find_elements(self.DRIVER.wd_by.CSS_SELECTOR, "a")

    _selected_link = None
    for _a in _anchors:
        _link = _a.get_attribute('href')
        if "/RNAfold/" in _link:
            _selected_link = _link

    self.config.rnafold.rna_submitted[_id] = (self.time_stamp(), _selected_link)
    self.update_config()

  def download_2d(self, _id, _url):
    if _id in self.config.rnafold.rna_collected.keys() or self.exists(f"{self.storage_path}/{_id}--result.html"):
        return

    _url_data = self.DRIVER.get_url(_url)

    _anchors = self.CHROME_INSTANCE.find_elements(self.DRIVER.wd_by.CSS_SELECTOR, "a")

    with open(f"{self.storage_path}/{_id}--result.html", "w+") as _fh:
        _fh.write(_url_data)

    for _a in _anchors:
        "Adobe SVG plugin"
        _link = _a.get_attribute('href')
        _fn = self.filename(_link, with_ext=True)
        if not _fn.startswith(_id):
            _fn = f"{_id}--{_fn}"
        _fp = f"{self.storage_path}/{_fn}"
        if _link.endswith(".eps") or _link.endswith(".fa") or _link.endswith(".out") or _link.endswith(".vienna") or _link.endswith(".ct"):
            self.URLREQ.urlretrieve(_link, _fp)

    self.config.rnafold.rna_collected[_id] = (self.time_stamp())
    self.update_config()

  def parse_2D(self, _rdata):
    self.re_energy = self.re_compile(r'-?\d+\.?\d+ kcal/mol')
    _ss_score_coll = []

    for _rna_id, _rna_seq in self.ProgressBar(_rdata):
        _files = self.find_files(self.storage_path)
        _macc_files = self.filter(_files, _rna_id)
        _result_files = self.filter(_macc_files, "--result")
        _score_parsed = {
          "rna_id": _rna_id,
          "rna_score_mfe": "NA",
          "rna_ss_centroid": "NA",
          "rna_energies": "NA",
        }
        if len(_result_files) > 0:
          _html = self.read_html(_result_files[0])

          _text_content = _html.find(True, {"id": "contentmain"}).text
          _energies = "|".join(self.re_energy.findall(_text_content))

          _mfe_ss = _html.find(True, {"id": "CENTROID_structure_span"}).find("pre").text.strip("1 ")
          _cent_ss = _html.find(True, {"id": "MFE_structure_span"}).find("pre").text.strip("1 ")

          _score_parsed = {
            "rna_id": _rna_id,
            "rna_score_mfe": _mfe_ss,
            "rna_ss_centroid": _cent_ss,
            "rna_energies": _energies,
          }
        _ss_score_coll.append(_score_parsed)

    self.df_score = self.DF(_ss_score_coll)
    self.df_score.to_pickle(self.get_path("rnafold.df.gz"))
    self.update_config()

  def get_2D(self, _rdata):
    # List of Tuples like ([rna_id, rna_seq], ...)
    self.init_scrap()
    for _id, _seq in self.ProgressBar(_rdata):
      self.submit_2d(_id, _seq)
    self.finish_scrap()

    self.time_sleep(20)
    self.init_scrap()
    for _rna_id, (_dt, _url) in self.ProgressBar(self.config.rnafold.rna_submitted.items()):
        self.download_2d(_rna_id, _url)
    self.finish_scrap()
    self.parse_2D(_rdata)

class Seq2D_CentroidFold(Seq2D_RNA):
  def __init__(self, *args, **kwargs):
    super(Seq2D_CentroidFold, self).__init__(**kwargs)
    self.__defaults = {
      "subversion": "centroidfold",
      "project_dir": "centroidfold",
      "url_server": "http://rtools.cbrc.jp/centroidfold/"
    }
    self.update_attributes(self, kwargs, self.__defaults)
    self.storage_path = self.validate_dir(self.get_path(self.project_dir))

  def submit_download_2d(self, _id, _seq):
    if _id in self.config.centroidfold.rna_collected.keys() or self.exists(f"{self.storage_path}/{_id}--structure.txt"):
        return

    _seq_info = f""">{_id}
{_seq.strip()}
"""

    _url_data = self.DRIVER.get_url(self.url_server)

    _query = self.CHROME_INSTANCE.find_elements(self.DRIVER.wd_by.NAME, "query")[0]
    _query.click()
    _query.send_keys(_seq_info)

    _select = self.DRIVER.wd_ui.Select(self.CHROME_INSTANCE.find_elements(self.DRIVER.wd_by.NAME, "gamma_1")[0])
    _select.select_by_value('16')

    _submit = self.CHROME_INSTANCE.find_elements(self.DRIVER.wd_by.CSS_SELECTOR, "input[type='submit']")[0]
    _submit.click()

    _new_url = self.CHROME_INSTANCE.current_url

    _anchors = self.CHROME_INSTANCE.find_elements(self.DRIVER.wd_by.CSS_SELECTOR, "div > a")

    for _a in _anchors:
        _link = _a.get_attribute('href')
        _fn = self.filename(_link, with_ext=True)
        _fp = f"{self.storage_path}/{_id}--{_fn}"
        if _link.endswith(".png") or _link.endswith(".txt"):
            self.URLREQ.urlretrieve(_link, _fp)

    self.config.centroidfold.rna_collected[_id] = {
        self.time_stamp(), _new_url
    }

    self.update_config()

  def parse_2D(self, _rdata):
    self.storage_path = self.get_path(self.project_dir)
    self.reg_score = self.re_compile(r'(.*)\s\((.*)\)')

    _ss_score_coll = []

    for _idx, _rna_seq in self.ProgressBar(_rdata):
        _text_rna = self.read_text(f"{self.storage_path}/{_idx}--structure.txt")
        _rna_ss, _rna_score = self.reg_score.findall(_text_rna[-1])[0]
        _score_parsed = {
            "rna_acc": _idx,
            "rna_ss": _rna_ss,
            "rna_score": _rna_score,
        }
        _ss_score_coll.append(_score_parsed)

    self.df_score = self.DF(_ss_score_coll)
    self.df_score.to_pickle(self.get_path("centroidfold.df.gz"))
    self.update_config()

  def get_2D(self, _rdata):
    # List of Tuples like ([rna_id, rna_seq], ...)
    self.init_scrap()
    for _id, _seq in self.ProgressBar(_rdata):
      self.submit_download_2d(_id, _seq)
    self.finish_scrap()
    self.parse_2D(_rdata)


class Seq3D_RNAComposer(SeqScrap):
  def __init__(self, *args, **kwargs):
    super(Seq3D_RNAComposer, self).__init__(**kwargs)
    self.__defaults = {
      "subversion": "rnacomposer",
      "project_dir": "rnacomposer",
      "url_server": "https://rnacomposer.cs.put.poznan.pl/",
      "url_workspace": "https://rnacomposer.cs.put.poznan.pl/Account/MyWorkspace",
      "requires_login": True,
    }
    self.update_attributes(self, kwargs, self.__defaults)
    self.storage_path = self.validate_dir(self.get_path(self.project_dir))

  def submit_batch_sequence(self, _id, _fs):
    if _id in self.config.RNAComposer.collected_set.keys():
      return

    _url_data = self.DRIVER.get_url(self.url_server)

    if self.requires_login:
      input("require user intervention for login.")
      self.requires_login = False

    _query = self.CHROME_INSTANCE.find_elements(self.DRIVER.wd_by.NAME, "content")[0]
    _query.clear()
    _query.click()
    _query.send_keys(_fs)

    _select = self.DRIVER.wd_ui.Select(self.CHROME_INSTANCE.find_elements(self.DRIVER.wd_by.NAME, "models")[0])
    _select.select_by_visible_text('5')

    _submit = self.CHROME_INSTANCE.find_elements(self.DRIVER.wd_by.NAME, "send")[0]
    _submit.click()

    self.config.RNAComposer.collected_set[_id] = {
      self.time_stamp()
    }

    self.update_config()

  def download_results(self, *args, **kwargs):
    # self.init_scrap(**kwargs)
    self.DRIVER.get_url(self.url_server)
    if self.requires_login:
      input("Check for login or other necessary actions.")

    self.DRIVER.get_url(self.url_workspace)
    _table_anchors = self.CHROME_INSTANCE.find_elements(self.DRIVER.wd_by.CSS_SELECTOR, "table.set_list_table a")

    # Remove all timeouts to stop refreshing the page
    self.CHROME_INSTANCE.execute_script("""var highestTimeoutId = setTimeout(";");
    for (var i = 0 ; i < highestTimeoutId ; i++) {
        clearTimeout(i);
    };""")

    _download_links = []

    self.validate_dir(self.get_path('Results--RNAComposer'))

    for _a in _table_anchors:
        _link = _a.get_attribute('href')
        if _a.text == 'Download':
          _download_links.append(_link)

    print(f"Download Links {len(_download_links)}/{len(_table_anchors)}", )
    self.clone_session()

    for _link in self.ProgressBar(_download_links):
        _set_id = _link.split("setId=")[-1]
        _file_path = self.get_path(f"Results--RNAComposer/{_set_id}.zip")

        if self.exists(_file_path):
            continue

        _res = self.Session_Clone.get(_link, stream=True)

        with open(_file_path, "wb") as _zf:
            _zf.write(_res.content)

        self.time_sleep(5)

    _table_rows = self.CHROME_INSTANCE.find_elements(self.DRIVER.wd_by.CSS_SELECTOR, "table.set_list_table tr")
    for _row in _table_rows:
        if "Download" in _row.text:
            _cb = _row.find_element(self.DRIVER.wd_by.CSS_SELECTOR, "input.id_checkbox")
            _set_id = _cb.get_attribute("value")
            _file_path = self.get_path(f"Results--RNAComposer/{_set_id}.zip")
            if self.exists(_file_path):
                _cb.click()

  def get_zip_file_info(self, *args, **kwargs):
    _processed_mirs = []
    _reg_mir = self.re_compile(r'\bMI\w+\d+')

    # Check all miR IDs in the ZIP files
    _zip_files = self.find_files(self.get_path(f"Results--RNAComposer"), "*.zip")

    for _zip in self.ProgressBar(_zip_files):
        self.path_zip = _zip
        _files = " ".join(self.list_zip_files(_zip, "names"))
        _files = _reg_mir.findall(_files)
        _processed_mirs.append({
            "zip_uid": self.filename(_zip),
            "n_files": len(_files),
            "mirs": list(set(_files)),
        })

    _zip_info = self.DF(_processed_mirs)
    return _zip_info

  def login(self, *args, **kwargs):
    self.DRIVER.get_url(self.url_server)

    _uname = args[0] if len(args) > 0 else kwargs.get("user")
    _pwd = args[1] if len(args) > 1 else kwargs.get("pwd")

    _form = self.CHROME_INSTANCE.find_element(self.DRIVER.wd_by.CSS_SELECTOR, "#menu form")

    _fld_usr = _form.find_element(self.DRIVER.wd_by.NAME, "j_username")
    _fld_usr.clear()
    _fld_usr.click()
    _fld_usr.send_keys(_uname)

    _fld_pwd = _form.find_element(self.DRIVER.wd_by.NAME, "j_password")
    _fld_pwd.clear()
    _fld_pwd.click()
    _fld_pwd.send_keys(_pwd)

    _xpath__submit = '//*[@id="menu"]/tbody/tr[2]/td/form/table/tbody/tr[3]/td/table/tbody/tr/td[1]/input'
    _fld_submit = _form.find_element(self.DRIVER.wd_by.XPATH, _xpath__submit)
    _fld_submit.click()


class Seq3D_3dRNA(SeqScrap):
  def __init__(self, *args, **kwargs):
    super(Seq3D_3dRNA, self).__init__(**kwargs)
    self.__defaults = {
      "subversion": "web3dRNA",
      "project_dir": "web3dRNA",
      "url_server": "http://biophy.hust.edu.cn/new",
      "ep_create": "%s/3dRNA/create",
      "ep_login": "%s/login",
      "ep_job_status": "%s/api/3dRNA/jobs/?page=%d",
      "ep_job_details": "%s/3dRNA/jobs",
      "ep_job_download": "%s/3dRNA/jobs/%s/download?name=%s",
      "requires_login": True,
      "path_3dRNA_status": "All-3dRNA-Result-Updates.df.gz",
      "result_types_3dRNA": ['all', 'par', 'log'],
    }
    self.__defaults.update(kwargs)
    self.update_attributes(self, self.__defaults)
    self.storage_path = self.validate_dir(self.get_path(self.project_dir))

  """
  seq_details = name = content
  num_models = name = models

  """

  def submit_sequence(self, _id, _seq, _ss):
    if _id in self.config.web3dRNA.collected_set.keys():
      return

    _url_data = self.DRIVER.get_url(self.ep_create % self.url_server)

    if self.requires_login:
      input("require user intervention for login.")
      self.requires_login = False

    _form = self.CHROME_INSTANCE.find_element(self.DRIVER.by_css, "form.el-form")

    _xpath__jobname = '//*[@id="pane-new_task"]/div/form/div[7]/div/div/div/input'
    _inp = _form.find_element(self.DRIVER.by_xpath, _xpath__jobname)
    _inp.clear()
    _inp.send_keys(_id)

    _xpath__seq = '//*[@id="pane-new_task"]/div/form/div[10]/div/div/div/textarea'
    _inp = _form.find_element(self.DRIVER.by_xpath, _xpath__seq)
    _inp.clear()
    _inp.send_keys(_seq)

    _xpath__ss = '//*[@id="pane-new_task"]/div/form/div[14]/div/div/div/textarea'
    _inp = _form.find_element(self.DRIVER.by_xpath, _xpath__ss)
    _inp.clear()
    _inp.send_keys(_ss)

    _xpath__sbmit_btn = '//*[@id="pane-new_task"]/div/form/div[22]/div/button[1]/span'
    _inp = _form.find_element(self.DRIVER.by_xpath, _xpath__sbmit_btn)
    _inp.click()

    self.config.web3dRNA.collected_set[_id] = {
      self.time_stamp()
    }

    self.update_config()

  def get_page_status(self, *args, **kwargs):
    _page_id = args[0] if len(args) > 0 else kwargs.get("page_id", 1)
    _page_url = self.ep_job_status % (self.url_server, _page_id)

    self.validate_dir(self.get_path("JSONS--3dRNA"))

    _json_file_path = self.get_path(f"JSONS--3dRNA/{_page_id}.pkl.gz")
    if self.exists(_json_file_path):
      return self.unpickle(_json_file_path)
    else:
      _res = self.Session_Clone.get(_page_url, stream=True)
      _json = _res.json()
      self.pickle(_json_file_path, _json)
      self.time_sleep(5)
      return _json

  def get_all_job_status(self, *args, **kwargs):
    if not self.exists(self.get_path(self.path_3dRNA_status)):
      return None
    else:
      return self.unpickle(self.get_path(self.path_3dRNA_status))

  def download_results(self, *args, **kwargs):
    self.clone_session()
    self.Session_Clone.headers.update({
        "X-XSRF-TOKEN": self.CHROME_INSTANCE.get_cookie('XSRF-TOKEN').get('value'), # Handle none
        })

    _json_first_page = self.get_page_status(1)
    _page_start = int(_json_first_page.get('from', 1))
    _page_last = int(_json_first_page.get('last_page', 1))

    self.validate_dir(self.get_path("Results--3dRNA"))

    for _pid in self.ProgressBar(range(_page_last, _page_start-1, -1), total=_page_last):
      _p_status = self.get_page_status(_pid)
      _df = self.DF(_p_status.get('data', []))

      for _, _job in _df.iterrows():
        if 'compressing' in _job.steps:
          self.download_result(_job.uuid)
      self.time_sleep(4)

  def download_result(self, *args, **kwargs):
    _uuid = args[0] if len(args) > 0 else kwargs.get("uuid")

    if not _uuid:
      return False

    if not self.Session_Clone:
      self.clone_session()

    for _res_typ in self.result_types_3dRNA:
      _res_url = self.ep_job_download % (self.url_server, _uuid, _res_typ)
      _res_path = self.get_path(f"Results--3dRNA/{_uuid}.{_res_typ}")

      if _res_path.endswith("all"):
        _res_path = f"{_res_path}.zip"

      if self.exists(_res_path):
        continue

      _res = self.Session_Clone.get(_res_url, stream=True)

      with open(_res_path, "wb") as _fp:
        _fp.write(_res.content)

    self.time_sleep(3)

  def login(self, *args, **kwargs):

    _uname = args[0] if len(args) > 0 else kwargs.get("user")
    _pwd = args[1] if len(args) > 1 else kwargs.get("pwd")

    self.DRIVER.get_url(self.ep_login % self.url_server)

    self.time_sleep(5)

    _fld_usr = self.CHROME_INSTANCE.find_element(self.DRIVER.by_id, "email")
    _fld_usr.clear()
    _fld_usr.click()
    _fld_usr.send_keys(_uname)

    _fld_pwd = self.CHROME_INSTANCE.find_element(self.DRIVER.by_id, "password")
    _fld_pwd.clear()
    _fld_pwd.click()
    _fld_pwd.send_keys(_pwd)

    _fld_pwd = self.CHROME_INSTANCE.find_element(self.DRIVER.by_id, "remember")
    _fld_pwd.click()

    _fld_captcha = self.CHROME_INSTANCE.find_element(self.DRIVER.by_id, "captcha")
    _fld_captcha.clear()

    _captcha = input("Enter Captcha: ")
    _fld_captcha.click()
    _fld_captcha.send_keys(_captcha)

    _xpath__submit = '//*[@id="content"]/div/div/div/div[2]/form/div[5]/div/button'
    _fld_submit = self.CHROME_INSTANCE.find_element(self.DRIVER.by_xpath, _xpath__submit)
    _fld_submit.click()

    self.requires_login = False
