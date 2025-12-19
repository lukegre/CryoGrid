# --- Load Environment Variables ---
# This looks for a file named .env in the same directory
-include .env
export

# --- Configuration ---
S3_ENDPOINT_URL := https://os.zhdk.cloud.switch.ch
S3_PATH_PREFIX  := s3://spi-pamir-cryogrid/cryogrid_runs-luke
RUNS_DIR        := /cluster/scratch/$(USER)/cryogrid-runs
REPO_URL        := https://github.com/lukegre/CryoGrid.git

# Set S3 integrity fixes
export AWS_REQUEST_CHECKSUM_CALCULATION := WHEN_REQUIRED
export AWS_RESPONSE_CHECKSUM_VALIDATION := WHEN_REQUIRED

# Logic for CRYOGRID_RUN_NAME (CLI 'name=' takes priority over .env/Env Var)
RUN_NAME   := $(or $(name),$(CRYOGRID_RUN_NAME))
LOCAL_PATH := $(RUNS_DIR)/$(RUN_NAME)
S3_PATH    := $(S3_PATH_PREFIX)/$(RUN_NAME)/

.PHONY: help init download upload submit check-env check-name

help: ## Show this help message
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-15s\033[0m %s\n", $$1, $$2}'

init: ## Clone repo and setup scratch symlinks
	@echo "Cloning CryoGrid repository..."
	git clone --depth 1 $(REPO_URL) $(HOME)/CryoGrid
	@echo "Setting up scratch directories..."
	mkdir -p $(RUNS_DIR)
	ln -snf $(RUNS_DIR) $(HOME)/cryogrid-runs

download: check-env check-name ## Download run files from S3
	@echo "Downloading $(RUN_NAME) from S3..."
	mkdir -p $(LOCAL_PATH)
	aws s3 sync $(S3_PATH) $(LOCAL_PATH) --endpoint-url $(S3_ENDPOINT_URL)

upload: check-env check-name ## Upload run results to S3
	@echo "Uploading $(RUN_NAME) to S3..."
	aws s3 sync $(LOCAL_PATH) $(S3_PATH) --endpoint-url $(S3_ENDPOINT_URL)

submit: check-name ## Submit SLURM job
	@echo "Submitting job for $(RUN_NAME)..."
	cd $(LOCAL_PATH) && sbatch sbatch_submit.sh

# --- Validation Helpers ---

check-env: ## Ensure AWS credentials are loaded
ifndef AWS_ACCESS_KEY_ID
	$(error AWS_ACCESS_KEY_ID is not set. Check your .env file)
endif
ifndef AWS_SECRET_ACCESS_KEY
	$(error AWS_SECRET_ACCESS_KEY is not set. Check your .env file)
endif

check-name: ## Ensure run name is provided
ifeq ($(RUN_NAME),)
	$(error RUN_NAME is not set. Use 'name=your_run' or set CRYOGRID_RUN_NAME in .env)
endif
