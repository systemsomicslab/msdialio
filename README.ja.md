# msdialio

## Rawデータの書き込み

`msdialio.raw` に `write_raw_arf`・`write_raw_arf2`・`write_raw_pai2`・
`write_raw_mddata`・`write_raw_dbs` を用意しています。
対応する `read_raw_*` の戻り値を渡してください。通常の `read_*` の戻り値は対象外です。

```python
from msdialio.raw import read_raw_pai2, write_raw_pai2

peaks = read_raw_pai2("sample.pai2")
# 必要な保存フィールドを編集
write_raw_pai2("sample_modified.pai2", peaks)
```

引数は `(出力先, rawデータ, *, overwrite=False)` です。パスとバイナリストリームに対応し、
呼び出し側のストリームは閉じません。パスは一時ファイル経由で保存し、既存ファイルへの
上書きには `overwrite=True` が必要です。ストリームは現在位置から書き込み、切り詰めず、
失敗時には途中まで書き込まれる場合があります。MessagePackの出力はメモリ上に保持します。

未知フィールド・数値キー・バイナリ・拡張型を保持します。float32は `float` のサブクラス
`msdialio.raw.Float32` として読み込み、float64と区別して保存します。float32の値を変更するときは
`Float32(新しい値)` を使ってください。演算や通常の `float` への変換では保存幅の情報が失われ、
通常の `float` はfloat64で保存されます。MS-DIALのfloat32読み込みはfloat64を受け付けません。
圧縮結果や整数の保存幅などのバイト単位の一致は保証しません。
C#のスキーマ、ファイル間のIDやパスの整合性は検証しません。
ARFは既存のブロック・オブジェクト構造・余白を保った値の編集を対象とし、行の追加・削除に
必要なLargeList構造の再構築は行いません。非標準の配列ヘッダ幅は保持しません。

DBSは `Storage` を正として `Entries["Storage"]` を再生成します。それ以外のエントリは
未知データも含めて元のバイト列を保存します。ZIPの日時などは保持しません。
Raw API内のエントリ名は `/` 区切りですが、書き込み時はMS-DIALの検索に合わせて
ZIP内部のエントリ名をWindows形式のバックスラッシュ区切りにします。
展開済みの参照レコードからDB本体を再構築する機能は含みません。

実ファイルの往復とバイナリ構造に加え、サンプルの編集結果をMS-DIALのGUIでも確認しています。
全バージョン・全データセットの互換性を保証するものではありません。
[C#互換性確認ツール](tools/msdial-interop/README.md)で実際の読み込み処理とLCMSデータセットの
初期化も確認できます。DCL・AEFの書き込みは対象外です。

[English](README.md) | 日本語

MS-DIALの `.arf`・`.arf2`・`.pai2`・`.dcl`・`.EIC.aef`・`.mddata`・`.dbs` を読むPythonパッケージです。
読み込みに加え、Raw APIからARF・ARF2・PAI2・MDDATA・DBSの書き込みに対応します。

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
DBSストリームにはseekが必要です。データはメモリ上に読み込みます。書き込みにはRaw APIを使用します。

### DBの選択と参照レコードの順次読み込み

`list_dbs` でDB一覧を確認してから、IDを指定して読み込めます。
ZIPのエントリ一覧と `Storage` だけを読み、DB本体や注釈器の実データは開きません。
返り値は `DataBaseID`・`Kind`・`DataBase`・`Pairs` を含む辞書のリストで、
メタデータの形式と並び順は `read_dbs` と同じです。
`References`・`Data`・化合物数は含めません。
DB本体の欠落や破損は一覧取得時には検証せず、本体を読み込む際に検出します。
パスまたはseek可能なバイナリストリームを受け取り、渡されたストリームは閉じません。

```python
from msdialio import list_dbs, read_dbs

path = "Dataset_Loaded.msp2.dbs"
catalog = list_dbs(path)
for database in catalog:
    print(database["DataBaseID"], database["Kind"], database["DataBase"])

if catalog:
    selected = read_dbs(path, database_ids=[catalog[0]["DataBaseID"]])
```

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
