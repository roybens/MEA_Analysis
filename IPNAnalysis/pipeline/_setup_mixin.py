# ==========================================================
# pipeline/_setup_mixin.py
# SetupMixin — logger, runtime controls, metadata parsing,
# and path-validation helpers used during MEAPipeline.__init__.
# ==========================================================

import logging
import os
import re
import sys
import configparser
from pathlib import Path


class SetupMixin:
    """Provides initialisation helpers for MEAPipeline."""

    # ------------------------------------------------------------------
    # Logger
    # ------------------------------------------------------------------

    def _setup_logger(self, log_file):
        logger = logging.getLogger(f"mea_{self.stream_id}")
        logger.setLevel(logging.DEBUG if self.verbose else logging.INFO)
        # Prevent duplicate handlers when running in a loop
        if not logger.handlers:
            formatter = logging.Formatter('[%(asctime)s] %(levelname)s: %(message)s')
            fh = logging.FileHandler(log_file, mode='a')
            fh.stream.write("\n" + "=" * 80 + "\n")
            fh.setFormatter(formatter)
            logger.addHandler(fh)
            ch = logging.StreamHandler(sys.stdout)
            ch.setFormatter(formatter)
            logger.addHandler(ch)
        return logger

    # ------------------------------------------------------------------
    # Runtime / resource controls
    # ------------------------------------------------------------------

    def _apply_runtime_controls(self):
        if self.cuda_visible_devices is not None:
            try:
                os.environ["CUDA_VISIBLE_DEVICES"] = str(self.cuda_visible_devices)
            except Exception:
                pass

    def _log_runtime_controls(self):
        def _env_or_none(name):
            value = os.environ.get(name)
            if value is None:
                return None
            token = str(value).strip()
            return token if token else None

        self.logger.info(
            "Runtime snapshot: pid=%s cpu_count=%s n_jobs=%s chunk_duration=%s",
            os.getpid(),
            os.cpu_count(),
            self.n_jobs,
            self.chunk_duration,
        )
        self.logger.info(
            "Runtime controls: cuda_visible_devices=%s",
            self.cuda_visible_devices,
        )
        self.logger.info(
            "Runtime env effective: CUDA_VISIBLE_DEVICES=%s",
            _env_or_none("CUDA_VISIBLE_DEVICES"),
        )

    # ------------------------------------------------------------------
    # Metadata
    # ------------------------------------------------------------------

    def _parse_metadata(self):
        meta = {
            'run_id': None,
            'chip_id': None,
            'project': None,
            'relative_pattern': (
                f"{self.file_path.parent.parent.name}"
                f"/{self.file_path.parent.name}"
                f"/{self.file_path.name}"
            ),
            'date': None,
            'well': None,
        }

        # Strategy A: Regex on path (fallback)
        try:
            path_str = str(self.file_path)
            match = re.search(r"/(\d+)/data.raw.h5", path_str)
            if match:
                meta['run_id'] = match.group(1)
            parts = path_str.split(os.sep)
            if len(parts) > 5:
                meta['relative_pattern'] = os.path.join(*parts[-6:-1])
                meta['project'] = parts[-6]
                meta['date'] = parts[-5]
                meta['chip_id'] = parts[-4]
                meta['well'] = self.stream_id
        except Exception:
            pass

        # Strategy B: .metadata file (overrides regex)
        meta_file = self.file_path.parent / ".metadata"
        if meta_file.exists():
            try:
                cfg = configparser.ConfigParser()
                cfg.read(meta_file, encoding='utf-8')
                if 'properties' in cfg:
                    meta['run_id'] = cfg['properties'].get('runid', meta.get('run_id'))
                    meta['project'] = cfg['properties'].get('project_title', meta.get('project'))
                if 'runtime' in cfg:
                    meta['chip_id'] = cfg['runtime'].get('chipid', meta.get('chip_id'))
            except Exception:
                pass

        return meta

    # ------------------------------------------------------------------
    # Path validation
    # ------------------------------------------------------------------

    def _validate_output_subdir_after_well(self, value):
        if value is None:
            return None

        token = str(value).strip()
        if not token:
            return None

        if "/" in token or "\\" in token:
            raise ValueError(
                "output_subdir_after_well must be a single directory name, not a path"
            )

        candidate = Path(token)
        if candidate.is_absolute() or token in (".", ".."):
            raise ValueError(
                "output_subdir_after_well must be a relative single directory name"
            )

        return token
