# ==========================================================
# pipeline/_checkpoint_mixin.py
# CheckpointMixin — checkpoint load, save, and skip logic.
# ==========================================================

import json
from datetime import datetime

from .stages import ProcessingStage, CHECKPOINT_SCHEMA_VERSION


class CheckpointMixin:
    """Manages pipeline checkpoints so interrupted runs can be resumed."""

    def _load_checkpoint(self):
        if self.checkpoint_file.exists() and not self.force_restart:
            with open(self.checkpoint_file, 'r') as f:
                state = json.load(f)

            # Migrate legacy stage numbering (schema v1 had no merge stage).
            try:
                schema_version = int(state.get("checkpoint_schema_version", 1))
            except Exception:
                schema_version = 1

            if schema_version < CHECKPOINT_SCHEMA_VERSION:
                try:
                    old_stage = int(state.get("stage", ProcessingStage.NOT_STARTED.value))
                except Exception:
                    old_stage = ProcessingStage.NOT_STARTED.value
                # v1 mapping: ANALYZER/ANALYZER_COMPLETE/REPORTS/REPORTS_COMPLETE
                # were 5/6/7/8; shift by +2 to account for new MERGE stages.
                if old_stage >= 5:
                    state["stage"] = old_stage + 2
                state["checkpoint_schema_version"] = CHECKPOINT_SCHEMA_VERSION

            return state

        return {
            'stage': ProcessingStage.NOT_STARTED.value,
            'checkpoint_schema_version': CHECKPOINT_SCHEMA_VERSION,
            'failed_stage': None,
            'last_updated': None,
            'run_id': self.run_id,
            'chip_id': self.chip_id,
            'well': self.well,
            'project': self.project_name,
            'date': self.date,
            'output_dir': str(self.output_dir),
            'error': None,
        }

    def _save_checkpoint(self, stage, **kwargs):
        self.state['stage'] = stage.value
        self.state['checkpoint_schema_version'] = CHECKPOINT_SCHEMA_VERSION
        self.state['last_updated'] = str(datetime.now())
        self.state.update(kwargs)
        with open(self.checkpoint_file, 'w') as f:
            json.dump(self.state, f, indent=2)
        self.logger.info(f"Checkpoint Saved: {stage.name}")

    def should_skip(self):
        if (
            self.state['stage'] == ProcessingStage.REPORTS_COMPLETE.value
            and not self.force_restart
        ):
            self.logger.info("Pipeline already completed. Skipping.")
            return True
        return False
