from setuptools import setup, find_packages

setup(
    name="MEA_Analysis",
    version="0.0.3",
    # description="Description of your project",
    # author="Your Name",
    # author_email="",
    # url="https://example.com/your-project",  # Replace with your project URL
    packages=find_packages(include=[
                            'AxonReconAnalysis',
                            'GUI',
                            'IPNAnalysis',
                            'MaxwellBiosystemsDeviceInterface',
                            'MEAProcessingLibrary',
                            'NetworkAnalysis',
                            'NeuronClassification',
                            'Organoid',
                            'WildtypeSegregation'
                            ]),
    install_requires=[
        # List dependencies here, e.g.,
        # 'numpy>=1.21.0',
    ],
)