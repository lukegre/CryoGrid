# --- Configuration ---
S3_ENDPOINT_URL := https://os.zhdk.cloud.switch.ch
S3_PATH_PREFIX  := s3://spi-pamir-cryogrid/cryogrid_runs-luke
RUNS_DIR        := /cluster/scratch/$(USER)/cryogrid-runs
REPO_URL        := https://github.com/lukegre/CryoGrid.git

# Set these in your shell or .env file
export AWS_REQUEST_CHECKSUM_CALCULATION := WHEN_REQUIRED
export AWS_RESPONSE_CHECKSUM_VALIDATION := WHEN_REQUIRED
export AWS_ACCESS_KEY_ID ?= your_key_here
export AWS_SECRET_ACCESS_KEY ?= your_secret_here

# Logic for CRYOGRID_RUN_NAME (Argument takes priority over Env Var)
# Usage: make upload name=my-run-name
RUN_NAME := $(or $(name),$(CRYOGRID_RUN_NAME))
LOCAL_PATH := $(RUNS_DIR)/$(RUN_NAME)
S3_PATH    := $(S3_PATH_PREFIX)/$(RUN_NAME)/

.PHONY: help init download upload submit check-name

help: ## Show this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'

init: ## Clone repo and setup scratch symlinks
	@echo "Cloning CryoGrid repository..."
	git clone --depth 1 $(REPO_URL) $(HOME)/CryoGrid
	@echo "Setting up scratch directories..."
	mkdir -p $(RUNS_DIR)
	ln -snf $(RUNS_DIR) $(HOME)/cryogrid-runs

download: check-name ## Download run files from S3 (Usage: make download name=RUN_ID)
	@echo "Downloading $(RUN_NAME) from S3..."
	mkdir -p $(LOCAL_PATH)
	aws s3 sync $(S3_PATH) $(LOCAL_PATH) --endpoint-url $(S3_ENDPOINT_URL)

upload: check-name ## Upload run results to S3 (Usage: make upload name=RUN_ID)
	@echo "Uploading $(RUN_NAME) to S3..."
	aws s3 sync $(LOCAL_PATH) $(S3_PATH) --endpoint-url $(S3_ENDPOINT_URL)

submit: check-name ## Submit SLURM job (Usage: make submit name=RUN_ID)
	@echo "Submitting job for $(RUN_NAME)..."
	cd $(LOCAL_PATH) && sbatch sbatch_submit.sh

check-name: ## Internal check to ensure name is provided
ifeq ($(RUN_NAME),)
	$(error RUN_NAME is not set. Use 'make <cmd> name=your_run' or set CRYOGRID_RUN_NAME env var)
endif