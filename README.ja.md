# MS-DIAL Parsers

[English](README.md) | [日本語](README.ja.md)

MS-DIAL が出力する `.arf`、`.arf2`、`.pai2`、`.dcl`、`.EIC.aef` を読む Python パッケージです。

## インストール

### リポジトリを clone して使う

NumPy を使わない場合:

```bash
python -m pip install .
```

uv を使う場合は、リポジトリ内で同期します。

```bash
uv sync
```

`.arf` の読み込み結果を数値行列に変換する utility `intensity_matrix` を使う場合だけ、追加依存 `matrix` を指定してください。NumPy なしでも他の読み込み API は使えます。関数を呼び出して NumPy がない場合は、インストール方法がエラーに表示されます。

pip:

```bash
python -m pip install ".[matrix]"
```

uv:

```bash
uv sync --extra matrix
```

### GitHub から直接インストールする

以下の `OWNER/REPO` は実際の GitHub リポジトリに置き換えてください。

通常インストール:

```bash
python -m pip install "msdial-parsers @ git+https://github.com/OWNER/REPO.git"
```

`intensity_matrix` を使う場合:

```bash
python -m pip install "msdial-parsers[matrix] @ git+https://github.com/OWNER/REPO.git"
```

uv（既存の仮想環境にインストールする場合）:

```bash
uv pip install "msdial-parsers @ git+https://github.com/OWNER/REPO.git"
```

`intensity_matrix` を使う場合:

```bash
uv pip install "msdial-parsers[matrix] @ git+https://github.com/OWNER/REPO.git"
```

## Python API

```python
from msdial_parsers import read_arf, peak_table

spots = read_arf("AlignmentResult_PeakProperties.arf")
rows = peak_table(spots)
```

対応する API は `read_arf`、`read_arf2`、`read_pai2`、`read_dcl`、`read_aef`、および各リーダーモジュールの補助関数です。

## コマンドライン

```bash
python -m msdial_parsers arf input.arf --json peaks.json
python -m msdial_parsers aef input.EIC.aef --csv spots.csv
```

形式は `aef`, `arf`, `arf2`, `dcl`, `pai2` から選びます。CLI の詳細は `python -m msdial_parsers --help` を参照してください。

## 開発

```bash
python -m pip install -e '.[test,matrix]'
python -m pytest
```

## ライセンス

MIT
