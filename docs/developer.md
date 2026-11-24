# Developer Guide

## Repository Layout

| Path | Purpose |
| --- | --- |
| `main.nf` | Entry point and workflow router. |
| `workflows/` | Top-level workflow implementations. |
| `subworkflows/` | Reusable workflow chains. |
| `modules/` | Process modules with tool calls and environments. |
| `bin/` | Python, Perl, and shell helper scripts. |
| `lib/` | Groovy helper classes for validation and pinned FooDMe2 asset access. |
| `conf/` | Runtime, resource, module, and reference configuration. |
| `docs/` | User and developer documentation. |
| `tests/` | Python unit tests. |
| `test/` | Nextflow module smoke-test harnesses and small fixtures. |

## Parameter Validation

`lib/WorkflowPipeline.groovy` owns high-level user validation:

- `--list_dbs`
- `--list_primers`
- mutually exclusive primer input modes
- required bounds for primer FASTA input (single file or directory)
- `--run_name` requirement when `--dbs` is used

Update this file when adding a user-facing parameter that changes valid run modes.

## Config Files

- `nextflow.config`: top-level defaults, plugin setup, profiles, includes.
- `conf/resources.config`: known reference databases and default installed paths.
- `conf/reference_sources.config`: source URLs for installable references.
- `conf/modules.config`: publish directories and process-specific arguments.
- `conf/base.config`: default resources and retry policy.
- `conf/modules/installation.config`: install workflow publish behavior.

## Adding A Module

1. Create `modules/<tool>/<task>/main.nf` or `modules/<tool>/main.nf`.
2. Add an `environment.yml` or container declaration.
3. Emit `versions.yml` and mix it into the parent workflow's version channel.
4. Add publish rules in `conf/modules.config` if outputs are user-facing.
5. Add a small test harness under `test/<module>/main.nf`.
6. Add Python unit tests when the module wraps logic in `bin/`.

`PrimerCatalog.groovy` keeps remote primer assets pinned and reproducible.
Analysis logic belongs in BarBeQuE modules; upstream workflow processes are not
imported.

## Testing

Python utility tests:

```bash
python3 -m unittest discover -s tests -p 'test*.py' -v
```

Full smoke-test runner:

```bash
bash run_all_tests.sh
```

The full runner needs Nextflow plus the selected software backend. In a minimal development shell, the Python suite may pass while module tests fail because Conda or containers are unavailable.

## Documentation

When changing workflow behavior, update:

- `docs/usage.md` for user-facing parameters
- `docs/pipeline.md` for top-level routing changes
- `docs/barbeque.md` for normal analysis changes
- `docs/output.md` for publish directory or output schema changes
- `nextflow_schema.json` for parameter help

Keep examples executable and avoid documenting parameters that are not wired into `main.nf` or `nextflow_schema.json`.
