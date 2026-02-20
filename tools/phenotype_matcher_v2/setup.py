"""
Setup script for the Phenotype Matcher v2 tool.
"""

from setuptools import setup, find_packages

setup(
    name="phenotype-matcher-v2",
    version="2.0.0",
    author="PAVS Team",
    description="Phenotype/disease matcher v2: formal algorithm with HPO, MONDO, OMIM, Orphanet",
    url="https://github.com/your-org/pavs",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "pronto>=2.5.0",
        "sentence-transformers>=2.2.0",
        "scikit-learn>=1.0.0",
        "requests>=2.28.0",
        "spacy>=3.0.0",
        "torch>=1.13.0",
        "rapidfuzz>=3.0.0",
        "pyahocorasick>=2.0.0",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0.0",
            "pytest-cov>=4.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "phenotype-matcher=phenotype_matcher_v2.cli:main",
        ],
    },
)
