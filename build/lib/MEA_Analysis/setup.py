# from setuptools import setup, find_packages

# setup(
#     name="MEA_Analysis",
#     version="0.0.3",
#     # description="Description of your project",
#     # author="Your Name",
#     # author_email="",
#     # url="https://example.com/your-project",  # Replace with your project URL
#     packages=find_packages(include=[
#                             'AxonReconAnalysis',
#                             'GUI',
#                             'IPNAnalysis',
#                             'MaxwellBiosystemsDeviceInterface',
#                             'MEAProcessingLibrary',
#                             'NetworkAnalysis',
#                             'NeuronClassification',
#                             'Organoid',
#                             'WildtypeSegregation'
#                             ]),
#     install_requires=[
#         # List dependencies here, e.g.,
#         # 'numpy>=1.21.0',
#     ],
# )

from setuptools import setup, find_packages
packages = find_packages(where="..") # start from parent of root so that the root is included
setup(
    name="MEA_Analysis",
    version="0.0.1",
    packages=packages,  # Explicitly search from the root
    #package_dir={"": "MEA_Analysis"},  # Map the root directory as the base
    #packages=find_packages(where=".."),  # Explicitly search from the root
    package_dir={"": ".."},  # Map the root directory as the base
    include_package_data=True,  # Ensure all package data is included
    install_requires=[],  # Add any dependencies if needed
)