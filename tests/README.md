# Repository Tests

These tests validate the repository structure and content without running any analysis.

## Running Tests

```bash
# From repository root
Rscript tests/parse_all_scripts.R
Rscript tests/expected_file_checks.R
Rscript tests/config_load_test.R
Rscript tests/dry_run_pipeline.R
```

Or run all tests:
```bash
for t in tests/*.R; do Rscript "$t"; done
```

## Test Descriptions

| Test | Purpose | Exit 0 |
|------|---------|--------|
| `parse_all_scripts.R` | Validates all R scripts parse without syntax errors | All parse OK |
| `expected_file_checks.R` | Checks required files exist (README, LICENSE, configs, etc.) | All present |
| `config_load_test.R` | Validates analysis_parameters.yml and paths.example.yml structure | All keys present |
| `dry_run_pipeline.R` | Lists all analysis steps and verifies scripts exist | All scripts found |

## Notes

- These tests do **not** execute any analysis
- They do **not** require raw data to be downloaded
- They validate repository structure only
