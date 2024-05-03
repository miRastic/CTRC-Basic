!pip install UtilityLib

## Sample Sheet to DataFrame
```py
from UtilityLib import ProjectManager
T0104 = ProjectManager(path_base=r"C:\Users\Abhimanyu\OneDrive\Desktop\GenomicData\TCGA")


_sample_sheet = T0104.pd_tsv(T0104.find_files(T0104.get_path(), "gdc_manifest*")[0])
_sample_sheet.columns = [_col.replace(' ', '_') for _col in _sample_sheet.columns]
T0104.pd_excel(T0104.get_path("sample-data.xlsx"), _sample_sheet, "Samples")
_sample_sheet.head(2)

```

## Metadata to DataFrame
```py
from UtilityLib import PM as T0104
T0104.set_project_paths([r"PATH_TO_DIR"]) # r"D:\datadrive\folder1"

_metadata_col_map = {
    "DataFormat": 'data_format',
    "FileName": 'file_name',
    "DataCategory": 'data_category',
    "DataType": 'data_type',
    "Experiment": 'experimental_strategy',
    "WorkflowType": 'workflow_type',
    "FileId": 'file_id',
    "Notes": 'annotations|0|notes',
    "annotations_case_id": 'annotations|0|case_id',
    "annotations_category": 'annotations|0|category',
    "annotations_status": 'annotations|0|status',
}

_metadata = T0104.json_to_df(T0104.find_files(T0104.get_path(), "metadata.cart.*.json")[0], _metadata_col_map)
T0104.pd_excel(T0104.get_path("sample-data.xlsx"), _sample_sheet, "MetaData")
_metadata.head(2)
```

## Expression File to DataFrame

```python

# For star count
if not isinstance(T0104.config.exp.star_files, (list, tuple)):
    T0104.config.exp.star_files = []

_expression_path = T0104.get_path('expression-star-count.pkl.gz')
if T0104.exists(_expression_path):
    _expression_df = T0104.unpickle(_expression_path)
else:
  T0104.validate_dir(T0104.get_path('Read-Chunks'))
  _chunk_size = 50
  _expression_files = T0104.find_files(T0104.get_path('downloads'), "*/*tsv")
  _chunk_list = list(enumerate(T0104.chunks(_expression_files, _chunk_size)))
  _chunk_paths = []
  for _chunk_idx, _file_set in T0104.PB(_chunk_list, total=len(_chunk_list)):
    _chunk_path = T0104.get_path(f'Read-Chunks/File-Chunk-{_chunk_idx}.pkl.gz')
    _chunk_paths.append(_chunk_path)
    if T0104.exists(_chunk_path):
        continue
    _expression_df = None
    for _f in _file_set:
        _fn = T0104.filename(_f, with_ext=True)
        if _f in T0104.config.exp.star_files:
            continue
        _df = T0104.pd_tsv(_f)
        _df['FileName'] = _fn
        _expression_df = _df if _expression_df is None else T0104.PD.concat([_expression_df, _df], axis=0)
        T0104.config.exp.star_files.append(_fn)

    _expression_df.to_pickle(_chunk_path)
    T0104.update_config()

  # Combine all the chunks
  _combined_exp_df = None
  for _cp in T0104.PB(_chunk_paths):
    _df = T0104.unpickle(_cp)
    _combined_exp_df = _df if _combined_exp_df is None else T0104.PD.concat([_combined_exp_df, _df], axis=0)
  _combined_exp_df.to_pickle(_expression_path)

```
