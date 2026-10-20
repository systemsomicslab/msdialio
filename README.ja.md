# msdialio

[English](README.md) | 日本語

MS-DIALの `.arf`・`.arf2`・`.pai2`・`.dcl`・`.EIC.aef`・`.mddata`・`.dbs` を読むPythonパッケージです。
現在は読み込みに対応し、MS-DIAL形式への書き込みは未実装です。

## インストール

リポジトリのルートで実行します。

```bash
python -m pip install .
# または
uv sync
```

## 通常のAPI

```python
from msdialio import read_arf, read_arf2, read_pai2, read_dcl, read_aef

spots = read_arf("alignment_PeakProperties.arf")
sample_peak = spots[0]["SamplePeaks"][0]
peaks = read_pai2("sample.pai2", slim=False)
```

- `read_arf2`: 保存されたスポット代表値の辞書のリスト。
- `read_arf`: スポット順の `SpotIndex` と、サンプル別辞書の `SamplePeaks`。
- `read_pai2`: 単一測定のピーク辞書のリスト。
- `read_dcl`: スペクトル辞書のリスト。`dcl_indices` で対象を指定できます。
- `read_aef`: サンプル別クロマトグラムを含むスポット辞書のリスト。`spot_ids`・`file_ids` で対象を指定できます。

ARF・ARF2・PAI2は既定の `slim=True` で主要項目を返し、`slim=False` で追加項目を含めます。
フィールド名はPascalCaseです。サンプルピークでは `RT`・`RI`・`DriftTime`・`Mz`・
`PeakHeightTop`・`PeakAreaAboveZero` を使用します。ARF・PAI2の未設定時間は `None` です。
ARF2は `MassCenter`・`HeightAverage` などのスポット固有名を保ち、時間座標は現状 `RT` のみです。
DCL・AEFのキーは現状snake_caseです。

ARFの外側に代表値やAlignment IDは生成しません。`SpotIndex` はファイル内の順番です。
代表値はARF2から取得し、ファイル間の対応は利用するデータで確認してください。
MS/MS実スペクトルは通常、対応するDCLファイルにあります。

各関数はパスかバイナリストリームを受け取り、渡されたストリームを閉じません。
DCL・AEFのストリームにはseekが必要です。

## Raw API

```python
from msdialio.raw import read_raw_arf, read_raw_arf2, read_raw_pai2
```

Raw関数はトップレベルには公開しません。圧縮を展開してMessagePackをデコードし、
入れ子・未知フィールド・数値の辞書キー・欠損値の元の表現を保持します。
フィールド名への変換やピーク行の抽出は行わず、`slim` 引数もありません。

- `read_raw_arf`: ブロックのリスト。各ブロックはデコードしたオブジェクトのリスト。
- `read_raw_arf2`: デコードしたオブジェクトのリスト。1個の場合も外側のリストを保持。
- `read_raw_pai2`: 単一のデコード済みオブジェクト。通常はピーク配列。

MessagePack文字列はstr、バイナリはbytesになります。不正なUTF-8は
`surrogateescape` で保持し、同じエラーハンドラーで再エンコードできます。
不完全なMessagePackはエラーにします。全データをメモリに保持するAPIです。
DCL・AEFはMessagePack形式ではないため、このRaw APIの対象外です。

## CLIと開発

```bash
python -m msdialio arf input.arf --json peaks.json
python -m msdialio aef input.EIC.aef --csv spots.csv
python -m pip install -e ".[test]"
python -m pytest
```

CLIはJSON・CSVに出力できます。入れ子のCSVセルはJSONで表現します。
詳細な契約と制限は英語READMEも参照してください。ライセンスはMITです。

## MS-DIAL 5 のデータセット・参照DB

```python
from msdialio import read_mddata, read_dbs

dataset = read_mddata("Dataset.mddata")
files = dataset["AnalysisFiles"]
project = dataset["Parameter"]["ProjectParam"]
databases = read_dbs("Dataset_Loaded.msp2.dbs")
references = databases[0]["References"]
```

`read_mddata` は1つの辞書を返し、解析ファイル一覧、アライメントファイル一覧、
基本解析設定、埋め込み参照情報を取得できます。保存されたパスの書き換えや、関連ファイルの
自動読み込みは行いません。参照DBが別ファイルに保存され、埋め込みの `MspDB` が空の場合もあります。

`read_dbs` はDBごとの辞書のリストを返します。`DataBaseID`・`Kind`・`DataBase`・
`Pairs`・`References` を含みます。化合物の `References` には `ScanID`・`Name`・
`PrecursorMz`・`Formula`・`SMILES`・`InChIKey`・`Spectrum` などを格納します。
`Spectrum` は保存時の配列で、各ピークの先頭2要素がm/zと強度です。
`ScanID` はDB内の識別子なので、注釈との結合時には `DataBaseID` も保持してください。

対象は **MS-DIAL 5** です。既知のファイル情報と基本設定にはC#のフィールド名を付けます。
列挙型は整数、`Formula`・`AdductType`・`ChromXs`・注釈器情報などの入れ子は保存時の形式を
保持します。未知の項目は数値キーのまま残し、測定方式固有の設定には名前を推測して付けません。
Proteomics・EAD lipid DBはレコードを解釈せず、`Data` にバイト列を保存し
`References=None` とします。現在のZIP構造と旧来の入れ子ZIP、Windowsの区切り文字に対応します。
DBSストリームにはseekが必要です。データはメモリ上に読み込みます。書き込みは未対応です。

### DBの選択と参照レコードの順次読み込み

```python
from contextlib import closing
from msdialio import read_dbs, iter_dbs_references

databases = read_dbs("Dataset_Loaded.msp2.dbs", database_ids=["MspDB"])

with closing(iter_dbs_references("Dataset_Loaded.msp2.dbs", database_id="MspDB")) as references:
    for reference in references:
        print(reference["Name"], reference["PrecursorMz"])
```

`database_ids` には1つのID文字列、またはIDのリストなどを指定できます。
省略すると全DB、空のリストなら0件です。存在しないIDは `KeyError` になります。
返却順は保存順です。異なる種類に同じIDがある場合、`read_dbs` は一致したDBをすべて返します。
読み込むのは `Storage` と指定DBのエントリだけで、対象外のDBや注釈器の実データは開きません。
Raw APIの `read_raw_dbs` は従来どおり全エントリを読み込みます。

`iter_dbs_references` は1つの化合物DBを指定し、`read_dbs` と同じ参照辞書を1件ずつ返します。
全レコードをリストに保持せず、圧縮ブロックごとに展開します。ただし、DB全体が1ブロックなら
その全体の展開は必要なので、メモリ削減は限定的な場合があります。化合物単位の直接アクセスは行いません。
途中で止めると後続ブロックは読み込まず、その部分の破損も検証しません。
エラーは反復時に発生します。途中終了時は上記の `closing` などでリソースを解放してください。
呼び出し側が渡したストリームは閉じません。Proteomics・EADはこのiteratorの対象外です。

```bash
python -m msdialio dbs Dataset_Loaded.msp2.dbs --database-id MspDB --json selected.json
```

`--database-id` は繰り返して複数指定できます。

```python
from msdialio.raw import read_raw_mddata, read_raw_dbs

storage = read_raw_mddata("Dataset.mddata")
archive = read_raw_dbs("Dataset_Loaded.msp2.dbs")
# archive["Storage"]: デコード済みメタデータ
# archive["Entries"]: ZIPエントリ名 -> 元のバイト列（未知のエントリも保持）
```

Raw APIでは未知の項目も保持します。ZIP内のパス区切りは `/` に統一します。

```bash
python -m msdialio mddata Dataset.mddata --json dataset.json
python -m msdialio dbs Dataset_Loaded.msp2.dbs --json databases.json
```

mddataのCSV出力は1行で、入れ子の項目はJSONです。CLIのJSON出力ではバイト列や拡張型は
文字列になります。正確な保持が必要な場合はPythonのRaw APIを使用してください。
フォーマットの根拠となるMS-DIAL本体のソースへのリンクは英語READMEに記載しています。
