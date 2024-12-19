# Contribution Guidelines for Oligo Designer Toolsuite

Thank you for your interest in contributing to Oligo Designer Toolsuite! Here's how you can effectively contribute:

## General Guidelines

1. **Fork and Clone**: Fork this repository and clone your fork to your local machine. Always work in a new branch.
2. **Commit and Push**: Commit changes to your branch and push them to your fork.
3. **Pull Request**: Create a pull request from your branch to this repository.

## Contribution Specifics

### Implementing Modules

1. **Inheritance from Base Class**: When creating a new function, make sure you put it in a class that inherit from the relevant base class. The base class offers a consistent structure shared across all classes in the sub-module. For example, a class designed for property filtering should inherit from `PropertyFilterBase` and must implement the `apply` method.

3. **Unit Tests**: Write unit tests for your contributions in the `test/test_<module_name>.py` file. Ensure they pass before submitting.

### Implementing Pipelines

1. **Use Pre-defined Modules**: Use our existing modules when building pipelines. If additional steps are needed, add them to the library, ensuring they align with our design.
2. **Raising Issues**: If in doubt about any aspect, raise an issue on the repository, we would be very happy to help!

### Documentation
Please make sure all contributions are documented using reStructuredText (reST). This standard is chosen for its compatibility with Sphinx, enabling us to generate well-structured and consistent project documentation. Here's an example of reST documentation format:
```
def function(arg1, arg2):
    """
    A brief description of the function.

    :param arg1: Description of arg1.
    :type arg1: int
    :param arg2: Description of arg2.
    :type arg2: str
    :returns: Description of return value.
    :rtype: bool
    """
    pass
```

## Thank you!
By adhering to these guidelines, you help maintain the quality of Oligo Designer Toolsuite. We appreciate your collaboration and contributions!
