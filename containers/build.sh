#!/usr/bin/env bash
# containers/build.sh
#
# Build all Docker images and optionally push to a registry.
#
# Usage:
#   ./containers/build.sh                      # build locally only
#   ./containers/build.sh --push               # build and push to GHCR
#   ./containers/build.sh --registry ghcr.io/myorg  # override registry
#
# On HPC (no Docker):
#   apptainer pull scclone-python.sif docker://ghcr.io/YOUR_ORG/scclone-python:1.0.0

set -euo pipefail

REGISTRY="${REGISTRY:-ghcr.io/TODO_ORG}"  # TODO: set to your registry
VERSION="${VERSION:-1.0.0}"
PUSH=false

for arg in "$@"; do
  case $arg in
    --push)     PUSH=true ;;
    --registry=*) REGISTRY="${arg#*=}" ;;
  esac
done

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

echo "=== Building scclone-python:${VERSION} ==="
docker build \
  -f "${ROOT}/containers/Dockerfile.python" \
  -t "scclone-python:${VERSION}" \
  -t "${REGISTRY}/scclone-python:${VERSION}" \
  -t "${REGISTRY}/scclone-python:latest" \
  "${ROOT}"

if [ "$PUSH" = "true" ]; then
  echo "=== Pushing scclone-python ==="
  docker push "${REGISTRY}/scclone-python:${VERSION}"
  docker push "${REGISTRY}/scclone-python:latest"
fi

echo ""
echo "Build complete."
echo ""
echo "To convert to Singularity SIF for HPC use:"
echo "  apptainer pull scclone-python.sif docker://${REGISTRY}/scclone-python:${VERSION}"
echo ""
echo "Or build directly from Docker daemon (if on local machine with Docker):"
echo "  apptainer build scclone-python.sif docker-daemon://scclone-python:${VERSION}"
