# MS-DIAL Parsers

[English](README.md) | [日本語](README.ja.md)

MS-DIAL が出力する `.arf`、`.arf2`、`.pai2`、`.dcl`、`.EIC.aef` を読む Python パッケージです。

## インストール

### リポジトリを clone して使う

```bash
python -m pip install .
```

uv を使う場合は、リポジトリ内で同期します。

```bash
uv sync
```

### GitHub から直接インストールする

以下の `OWNER/REPO` は実際の GitHub リポジトリに置き換えてください。

通常インストール:

```bash
python -m pip install "msdial-parsers @ git+https://github.com/OWNER/REPO.git"
```

uv（既存の仮想環境にインストールする場合）:

```bash
uv pip install "msdial-parsers @ git+https://github.com/OWNER/REPO.git"
```

## Python API

```python
from msdial_parsers import read_arf

spots = read_arf("AlignmentResult_PeakProperties.arf")
```

パッケージのトップレベルAPIは読み込み関数です。`read_arf`、`read_arf2`、`read_pai2`、`read_dcl`、`read_aef` を利用できます。`read_aef` は `spot_ids` と `file_ids` を指定して、対象のクロマトグラムだけを直接読み込めます。

## コマンドライン

```bash
python -m msdial_parsers arf input.arf --json peaks.json
python -m msdial_parsers aef input.EIC.aef --csv spots.csv
```

形式は `aef`, `arf`, `arf2`, `dcl`, `pai2` から選びます。CLI の詳細は `python -m msdial_parsers --help` を参照してください。

## 開発

```bash
python -m pip install -e '.[test]'
python -m pytest
```

## ライセンス

MIT
