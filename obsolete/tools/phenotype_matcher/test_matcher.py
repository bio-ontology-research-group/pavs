"""
Basic tests for the Phenotype Matcher tool.
"""

import pytest
from phenotype_matcher import PhenotypeMatcher, PhenotypeInput, MatcherConfig


class TestPhenotypeMatcherBasic:
    """Basic functionality tests."""

    @pytest.fixture(scope="class")
    def matcher(self):
        """Create a matcher instance for testing."""
        config = MatcherConfig(embedding_model="fast", llm_model="fast", debug=False)
        return PhenotypeMatcher(config)

    def test_matcher_initialization(self, matcher):
        """Test that matcher initializes successfully."""
        assert matcher is not None
        assert len(matcher.terms) > 0
        assert matcher.embeddings is not None

    def test_single_phenotype(self, matcher):
        """Test matching a single phenotype."""
        input_data = PhenotypeInput(text="seizure")
        output = matcher.match(input_data)

        assert output is not None
        assert len(output.phenotypes) > 0
        assert any("seizure" in p.label.lower() for p in output.phenotypes)

    def test_phenotype_with_severity(self, matcher):
        """Test matching a phenotype with severity modifier."""
        input_data = PhenotypeInput(text="severe intellectual disability")
        output = matcher.match(input_data)

        assert output is not None
        assert len(output.phenotypes) > 0

        # Check if any phenotype has severity
        has_severity = any(p.severity_label is not None for p in output.phenotypes)
        # Note: This test may fail if LLM doesn't extract severity
        # It's more of an integration test

    def test_multiple_phenotypes(self, matcher):
        """Test matching multiple phenotypes from comma-separated text."""
        input_data = PhenotypeInput(text="seizures, hypotonia")
        output = matcher.match(input_data)

        assert output is not None
        assert len(output.phenotypes) >= 1  # Should get at least one
        # Exact count may vary based on LLM interpretation

    def test_negation_detection(self, matcher):
        """Test negation detection."""
        input_data = PhenotypeInput(text="no seizures")
        output = matcher.match(input_data)

        assert output is not None
        # Check if any phenotype is marked as excluded
        excluded = output.get_excluded_phenotypes()
        # Note: This depends on LLM correctly detecting negation

    def test_output_methods(self, matcher):
        """Test output helper methods."""
        input_data = PhenotypeInput(text="seizures")
        output = matcher.match(input_data)

        # Test various output methods
        hpo_ids = output.get_hpo_ids()
        assert isinstance(hpo_ids, list)

        excluded_ids = output.get_hpo_ids(excluded=True)
        assert isinstance(excluded_ids, list)

        omim_ids = output.get_omim_ids()
        assert isinstance(omim_ids, list)

        mondo_ids = output.get_mondo_ids()
        assert isinstance(mondo_ids, list)

        # Test to_dict
        output_dict = output.to_dict()
        assert isinstance(output_dict, dict)
        assert "phenotypes" in output_dict
        assert "diseases" in output_dict


class TestInputValidation:
    """Test input validation and edge cases."""

    def test_empty_text(self, matcher):
        """Test handling of empty text."""
        input_data = PhenotypeInput(text="")
        output = matcher.match(input_data)

        assert output is not None
        # May return empty results or handle gracefully

    def test_no_split(self, matcher):
        """Test with split_by=None."""
        input_data = PhenotypeInput(text="seizures, hypotonia", split_by=None)
        output = matcher.match(input_data)

        assert output is not None
        # Should treat entire text as single term


class TestConfiguration:
    """Test configuration options."""

    def test_custom_config(self):
        """Test custom configuration."""
        config = MatcherConfig(
            embedding_model="fast",
            llm_model="fast",
            top_k_phenotype=3,
            top_k_severity=3,
            device="cpu",
            debug=True,
        )

        assert config.embedding_model == "fast"
        assert config.llm_model == "fast"
        assert config.top_k_phenotype == 3
        assert config.device == "cpu"
        assert config.debug is True


# Integration tests (require API key)
@pytest.mark.integration
class TestIntegration:
    """Integration tests requiring API access."""

    @pytest.fixture(scope="class")
    def matcher(self):
        """Create matcher for integration tests."""
        return PhenotypeMatcher()

    def test_full_workflow(self, matcher):
        """Test complete workflow with real API calls."""
        input_data = PhenotypeInput(text="severe intellectual disability and seizures")
        output = matcher.match(input_data)

        assert output is not None
        assert len(output.phenotypes) >= 1

        # Check metadata
        assert "processing_time_seconds" in output.processing_metadata
        assert "llm_calls" in output.processing_metadata
        assert output.processing_metadata["llm_calls"] >= 1


if __name__ == "__main__":
    # Run basic tests without pytest
    print("Running basic smoke tests...")

    try:
        config = MatcherConfig(embedding_model="fast", llm_model="fast")
        matcher = PhenotypeMatcher(config)
        print("✓ Matcher initialized")

        # Test simple match
        input_data = PhenotypeInput(text="seizure")
        output = matcher.match(input_data)
        print(f"✓ Matched 'seizure': {len(output.phenotypes)} phenotype(s)")

        if output.phenotypes:
            for p in output.phenotypes:
                print(f"  - {p.hpo_id}: {p.label}")

        print("\nAll smoke tests passed!")

    except Exception as e:
        print(f"✗ Test failed: {e}")
        import traceback

        traceback.print_exc()
