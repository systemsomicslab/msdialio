# msdialio

[English](README.md) | 日本語

MS-DIALの `.arf`・`.arf2`・`.pai2`・`.dcl`・`.EIC.aef` を読むPythonパッケージです。
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
