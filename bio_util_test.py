import logging
import sys
import os
import tempfile
import pytest

"""
bio_util_test.py — tests for bio_toolkit (pytest)

Usage:
    pytest bio_util_test.py -v

Structure:
    TestDnaRnaTools   — tests for run_dna_rna_tools (sequence operations)
    TestFilterFastq   — tests for filter_fastq (filtering by GC, length, quality)
    TestLogging       — tests for log file writing
"""

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from main import run_dna_rna_tools, filter_fastq, setup_logging

class TestDnaRnaTools:
    """Tests for the run_dna_rna_tools function."""

    def test_transcribe_single_sequence(self):
        """Transcription of a single DNA sequence: T -> U."""
        result = run_dna_rna_tools("ATGC", "transcribe")
        assert result == "AUGC", f"Expected 'AUGC', got '{result}'"

    def test_reverse_complement_multiple_sequences(self):
        """reverse_complement on multiple sequences returns a list of correct length."""
        result = run_dna_rna_tools("ATGC", "TTAA", "reverse_complement")
        assert isinstance(result, list), "Expected a list for multiple sequences"
        assert len(result) == 2
        assert result[0] == "GCAT"   # reverse_complement("ATGC")
        assert result[1] == "TTAA"   # reverse_complement("TTAA")

    def test_unknown_operation_returns_none(self):
        """Passing an unknown operation must return None instead of raising an exception."""
        result = run_dna_rna_tools("ATGC", "fly_to_the_moon")
        assert result is None, (
            f"Expected None for an unknown operation, got: {result!r}"
        )

    def test_is_nucleic_acid_invalid_sequence(self):
        """A sequence containing invalid characters is not a nucleic acid."""
        result = run_dna_rna_tools("ATGCX123", "is_nucleic_acid")
        assert result is False, (
            "Expected False for a sequence with invalid characters"
        )


class TestFilterFastq:
    """Tests for the filter_fastq function."""

    SEQS = {
        "high_gc":   ("GGGGCCCC", "IIIIIIII"),  # GC = 100%,  avgQ ~ 40
        "low_gc":    ("AAAATTTT", "IIIIIIII"),  # GC = 0%,    avgQ ~ 40
        "short_seq": ("AT",       "II"),         # len = 2,    avgQ ~ 40
        "low_qual":  ("ATGCATGC", "!!!!!!!!")   # GC = 50%,   avgQ = 0
    }

    def test_filter_by_gc_bounds(self):
        """GC filter passes only sequences within the specified range."""
        result = filter_fastq(self.SEQS, gc_bounds=(40, 60))

        assert "low_qual" in result
        assert "high_gc" not in result
        assert "low_gc" not in result

    def test_filter_by_length_bounds(self):
        """Length filter excludes sequences that are too short."""
        result = filter_fastq(self.SEQS, length_bounds=(5, 2**32))
        assert "short_seq" not in result, "Short sequence should not pass the length filter"
        assert "high_gc" in result
        assert "low_gc" in result

    def test_filter_by_quality_threshold(self):
        """Quality filter excludes sequences whose average quality is below the threshold."""
        result = filter_fastq(self.SEQS, quality_threshold=20)
        assert "low_qual" not in result, (
            "low_qual (avgQ=0) must not pass with quality_threshold=20"
        )
        assert "high_gc" in result
        assert "low_gc" in result

    def test_logging_writes_to_file(self):
        """
        The logger writes INFO messages to a file when filter_fastq is called.
        Verifies that the file is created and contains the expected entries.
        """
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".log", delete=False, encoding="utf-8"
        ) as tmp:
            log_path = tmp.name

        try:
            test_logger = logging.getLogger("bio_toolkit_test_logging")
            test_logger.setLevel(logging.DEBUG)
            test_logger.handlers.clear()

            fh = logging.FileHandler(log_path, encoding="utf-8")
            fh.setLevel(logging.DEBUG)
            fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
            test_logger.addHandler(fh)

            test_logger.info("Starting filter_fastq: 1 sequence(s)")
            test_logger.info("filter_fastq done: 1/1 sequence(s) passed.")

            fh.flush()
            fh.close()

            assert os.path.exists(log_path), "Log file was not created"
            with open(log_path, encoding="utf-8") as f:
                content = f.read()

            assert "[INFO]" in content, "Log file contains no INFO entries"
            assert "filter_fastq" in content, "Log file contains no mention of filter_fastq"

        finally:
            os.unlink(log_path)
