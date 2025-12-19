# --- Load Environment Variables ---
-include .env
export

# --- Configuration ---
S3_ENDPOINT_URL := https://os.zhdk.cloud.switch.ch
S3_PATH_PREFIX  := s3://spi-pamir-cryogrid/cryogrid_runs-luke
RUNS_DIR        := /cluster/scratch/$(USER)/cryogrid-runs
REPO_URL        := https://github.com/lukegre/CryoGrid.git

# S3 integrity fixes for newer AWS CLI versions
export AWS_REQUEST_CHECKSUM_CALCULATION := WHEN_REQUIRED
export AWS_RESPONSE_CHECKSUM_VALIDATION := WHEN_REQUIRED

# Variable priority: 1. CLI (name=x) 2. .env (CRYOGRID_RUN_NAME)
RUN_NAME   := $(or $(name),$(CRYOGRID_RUN_NAME))
LOCAL_PATH := $(RUNS_DIR)/$(RUN_NAME)
S3_PATH    := $(S3_PATH_PREFIX)/$(RUN_NAME)/

.PHONY: help init download upload submit check-env check-name

help: ## Show this help message
	@echo "Usage: make [target] [name=run-name]"
	@echo ""
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

init: ## Setup repo and scratch symlinks
	@echo "Cloning CryoGrid..."
	git clone --depth 1 $(REPO_URL) $(HOME)/CryoGrid
	mkdir -p $(RUNS_DIR)
	ln -snf $(RUNS_DIR) $(HOME)/cryogrid-runs

download: check-env check-name ## Sync files from S3 to local scratch
	@echo "Downloading $(RUN_NAME)..."
	mkdir -p $(LOCAL_PATH)
	aws s3 sync $(S3_PATH) $(LOCAL_PATH) --endpoint-url $(S3_ENDPOINT_URL)

upload: check-env check-name ## Sync local scratch results to S3
	@echo "Uploading $(RUN_NAME)..."
	aws s3 sync $(LOCAL_PATH) $(S3_PATH) --endpoint-url $(S3_ENDPOINT_URL)

submit: check-name ## Submit the SLURM job
	@echo "Submitting $(RUN_NAME)..."
	cd $(LOCAL_PATH) && sbatch sbatch_submit.sh

# --- Guards ---

check-env:
ifndef AWS_ACCESS_KEY_ID
	$(error AWS_ACCESS_KEY_ID is missing in .env)
endif

check-name:
ifeq ($(RUN_NAME),)
	$(error RUN_NAME is not set. Use name=... or set CRYOGRID_RUN_NAME in .env)
endif