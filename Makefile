# --- Load Environment Variables ---
-include .env
export

# --- Configuration ---
# S3 paths
S3_ENDPOINT_URL := https://os.zhdk.cloud.switch.ch
S3_PATH_PREFIX  := s3://spi-pamir-cryogrid/cryogrid_runs-luke
S3_PATH_FORCING := s3://spi-pamir-cryogrid/forcing/era5-pamirs-1990_2024-v251113.mat
# local paths
RUNS_DIR        := /cluster/scratch/$(USER)/cryogrid-runs
FORCING_DIR     := ${HOME}/cryogrid-forcing/

# S3 integrity fixes
export AWS_REQUEST_CHECKSUM_CALCULATION := WHEN_REQUIRED
export AWS_RESPONSE_CHECKSUM_VALIDATION := WHEN_REQUIRED

# Variable priority: 1. CLI (name=x) 2. .env (CRYOGRID_RUN_NAME)
RUN_NAME   := $(or $(name),$(CRYOGRID_RUN_NAME))
LOCAL_PATH := $(RUNS_DIR)/$(RUN_NAME)
S3_PATH    := $(S3_PATH_PREFIX)/$(RUN_NAME)/

.PHONY: help init install-aws download upload submit check-env check-name check-aws

help: ## Show this help message
	@echo "Usage: make [target] [name=run-name]"
	@grep -hE '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

install-aws: ## Install AWS CLI v2 to ~/.local/bin
	@echo "Downloading AWS CLI installer..."
	curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
	unzip -q awscliv2.zip
	@echo "Installing to $(HOME)/.local/bin..."
	./aws/install -i $(HOME)/.aws-cli -b $(HOME)/.local/bin --update
	@rm -rf aws awscliv2.zip
	@echo "Installation complete. Ensure $(HOME)/.local/bin is in your PATH."

forcing: check-aws check-env check-name ## Force download by removing local copy first
	@echo "Download ERA5 forcing for Pamirs..."
	mkdir -p $(FORCING_DIR)
	aws s3 sync $(S3_PATH_FORCING) $(FORCING_DIR) --endpoint-url $(S3_ENDPOINT_URL)

init: ## Setup scratch symlinks
	mkdir -p $(RUNS_DIR)
	ln -snf $(RUNS_DIR) $(HOME)/cryogrid-runs
	${MAKE} install-aws
	${MAKE} download-forcing

download: check-aws check-env check-name ## Sync files from S3 to local scratch
	@echo "Downloading $(RUN_NAME)..."
	mkdir -p $(LOCAL_PATH)
	aws s3 sync $(S3_PATH) $(LOCAL_PATH) --endpoint-url $(S3_ENDPOINT_URL)

upload: check-aws check-env check-name ## Sync local scratch results to S3
	@echo "Uploading $(RUN_NAME)..."
	aws s3 sync $(LOCAL_PATH) $(S3_PATH) --endpoint-url $(S3_ENDPOINT_URL)

submit: check-name ## Submit the SLURM job
	@echo "Submitting $(RUN_NAME)..."
	cd $(LOCAL_PATH) && sbatch sbatch_submit.sh

# --- Guards & Helpers ---

check-aws:
	@command -v aws >/dev/null 2>&1 || { echo >&2 "Error: AWS CLI not found. Run 'make install-aws' first."; exit 1; }

check-env:
ifndef AWS_ACCESS_KEY_ID
	$(error AWS_ACCESS_KEY_ID is missing in .env)
endif

check-name:
ifeq ($(RUN_NAME),)
	$(error RUN_NAME is not set. Use name=... or set CRYOGRID_RUN_NAME in .env)
endif