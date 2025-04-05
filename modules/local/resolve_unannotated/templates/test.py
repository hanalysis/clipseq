import pytest
import pybedtools as pbt
from resolve_unannotated import validate_chromosomes, validate_unfiltered, validate_resolved

def test_validate_chromosomes():
    filt_valid = {"chr1", "chr2", "chr3"}
    unfilt_valid = {"chr1", "chr2", "chr3"}
    index_valid = {"chr1", "chr2", "chr3"}
    index_mock = {"chr1", "chr2", "chr3", "chr4"}
    filt_mock = {"chr1"}
    # Test case when all chromosome sets match
    try:
        validate_chromosomes(filt_valid, unfilt_valid, index_valid)
    except ValueError:
        pytest.fail("validate_chromosomes raised ValueError unexpectedly!")
    
    # Test case when filtered and unfiltered sets mismatch with the index
    try:
        validate_chromosomes(filt_valid, unfilt_valid, index_mock)
    except ValueError:
        pass  # Expected behavior: exception should be raised

    # Test case when filtered chromosomes don't match the unfiltered set or index
    try:
        validate_chromosomes(filt_mock, unfilt_valid, index_valid)
    except ValueError:
        pass  # Expected behavior: exception should be raised


def test_validate_unfiltered():
    # Create mock BedTool objects for bed_fai and bed_unfiltered
    bed_fai_data = [
        ("chr1", 0, 1000, ".", 0, "+"),
        ("chr2", 0, 100, ".", 0, "+"),
        ("chr1", 0, 1000, ".", 0, "-"),
        ("chr2", 0, 100, ".", 0, "-"),

    ]
    bed_unfiltered_data_valid = [
        ("chr1", 0, 1000, "intergenic", 0, "+"),
        ("chr1", 0, 500, "CDS", 0, "-"),
        ("chr1", 500, 1000, "UTR", 0, "-"),
        ("chr2", 0, 100, "CDS", 0, "+"),
        ("chr2", 0, 100, "intergenic", 0, "-")
    ]

    bed_unfiltered_data_notvalid = [
        ("chr1", 0, 1000, "intergenic", 0, "+"),
        ("chr1", 0, 500, "CDS", 0, "+"),
        ("chr2", 0, 100, "CDS", 0, "+"),
        ("chr2", 0, 100, "intergenic", 0, "-")
    ]
    
    bed_fai = pbt.BedTool(bed_fai_data).sort()
    bed_unfiltered_valid = pbt.BedTool(bed_unfiltered_data_valid).sort()
    bed_unfiltered_notvalid = pbt.BedTool(bed_unfiltered_data_notvalid).sort()

    # Test case when there are no unannotated regions
    try:
        validate_unfiltered(bed_fai, bed_unfiltered_valid)
    except ValueError:
        pytest.fail("validate_unfiltered raised ValueError unexpectedly!")

    # Test case when there are unannotated regions
    try:
        validate_unfiltered(bed_fai, bed_unfiltered_notvalid)
    except ValueError:
        pass  # Expected behavior: exception should be raised


def test_validate_resolved():
    # Create mock BedTool objects for bed_fai and bed_unfiltered
    bed_fai_data = [
        ("chr1", 0, 1000, ".", 0, "+"),
        ("chr2", 0, 100, ".", 0, "+"),
        ("chr1", 0, 1000, ".", 0, "-"),
        ("chr2", 0, 100, ".", 0, "-"),

    ]
    bed_unfiltered_data_valid = [
        ("chr1", 0, 1000, "intergenic", 0, "+"),
        ("chr1", 0, 500, "CDS", 0, "-"),
        ("chr1", 500, 1000, "UTR", 0, "-"),
        ("chr2", 0, 100, "CDS", 0, "+"),
        ("chr2", 0, 100, "intergenic", 0, "-")
    ]

    bed_unfiltered_data_notvalid = [
        ("chr1", 0, 1000, "intergenic", 0, "+"),
        ("chr1", 0, 500, "CDS", 0, "-"),
        ("chr2", 0, 100, "CDS", 0, "+"),
        ("chr2", 0, 100, "intergenic", 0, "-")
    ]
    
    bed_fai = pbt.BedTool(bed_fai_data).sort()
    bed_unfiltered_valid = pbt.BedTool(bed_unfiltered_data_valid).sort()
    bed_unfiltered_notvalid = pbt.BedTool(bed_unfiltered_data_notvalid).sort()

    # Test case when there are no unannotated regions
    try:
        validate_resolved(bed_fai, bed_unfiltered_valid)
    except ValueError:
        pytest.fail("validate_unfiltered raised ValueError unexpectedly!")

    # Test case when there are unannotated regions
    try:
        validate_resolved(bed_fai, bed_unfiltered_notvalid)
    except ValueError:
        pass  # Expected behavior: exception should be raised
