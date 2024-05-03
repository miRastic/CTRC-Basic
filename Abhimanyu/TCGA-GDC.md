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
_expression_files = T0104.find_files(T0104.get_path(), "*/*tsv")

_expression_path = T0104.get_path('expression.pkl.gz')
if T0104.exists(_expression_path):
    _expression_df = T0104.unpickle(_expression_path)
else:
    _expression_df = None
    for _f in T0104.PB(_expression_files):
        _fn = T0104.filename(_f, with_ext=True)
        _df = T0104.PD.read_csv(_f, sep="\t", header=1)
        _df = _df[~_df.gene_id.str.startswith('N_')].copy()

        _df['FileName'] = _fn
        _expression_df = _df if _expression_df is None else T0104.PD.concat([_expression_df, _df], axis=0)

# print(_expression_df.head(2))

```
