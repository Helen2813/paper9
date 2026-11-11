# P9-DARC v0.1 — инструкция по запуску

**Назначение:** начать разработку универсального метода Paper 9, не изменяя исследования v001–v071, результаты и рукопись.

## Что создано

1. `METHOD_SPEC_v0.1.md` — формальная спецификация, уравнения и ограничения.
2. `BENCHMARK_PLAN_v0.1.md` — план внешней калибровки и сравнений.
3. `NOVELTY_AND_BOUNDARIES.md` — отличие от CELLECT/scDRS/CATCH и границы Paper 7/8.
4. `p9darc/` — минимальный исполняемый исследовательский прототип с семью тестами.

## Установка на MacBook

1. Сохранить архив `P9_DARC_methods_v001.zip` в `~/Downloads`.
2. Открыть Terminal и вставить одним блоком:

```bash
set -e
cd "$HOME/Desktop/paper9"
unzip -n "$HOME/Downloads/P9_DARC_methods_v001.zip" -d .
PYTHONPATH="$PWD/methods" python3 -m unittest discover -s methods/p9darc/tests -v
if [ ! -e results/methods/p9darc_v001_demo ]; then
  PYTHONPATH="$PWD/methods" python3 -m p9darc.cli demo --out results/methods/p9darc_v001_demo
else
  echo '[SKIP] Demo output already exists; no overwrite.'
fi
git status --short
```

## Что должно быть видно

`Ran 7 tests ... OK` и строки `[P9-DARC] PROTOTYPE ONLY / NO GWAS INFERENCE`, `[STATUS] NOT_CALIBRATED`.

Отчёт: `~/Desktop/paper9/results/methods/p9darc_v001_demo/audit_report.json`.

**Важно:** это тест на искусственных данных. Средние MAGMA gene-z в демо не корректируют LD, не являются p-value и не могут использоваться для утверждения о причине заболевания. Никаких реальных генетических результатов эти команды не создают. Для S1 понадобится подключение официального MAGMA/S-LDSC и замороженных выходных файлов Paper 9 в режиме чтения.

## Безопасность репозитория

Архив раскладывается в `methods/`, ключ `unzip -n` запрещает перезаписывать совпадающие файлы. Не создаётся Git commit, не модифицируются `docs/plan/CHANGELOG.md`, `scripts/p9_v054*`, `scripts/p9_v071*`, `data/` или результаты прошлых запусков.

## Следующая версия после этого шага

`methods/S1_DATA_ADAPTER_SPEC` и оригинальные MAGMA/S-LDSC wrappers по реальным фиксированным файлам, затем **независимый** методологический benchmark. Не выдавать заранее `STATE_SUPPORTED` до отдельной калибровки.
