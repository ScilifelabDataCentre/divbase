# Contributing to DivBase

This guide describes how to add development contributions to DivBase, which we welcome from anyone.

If you would like to suggest a more substantial change or new feature, please [create a GitHub Issue](https://github.com/ScilifelabDataCentre/divbase/issues) first so we can discuss it.

## Ways you can contribute

- **Bug Reports**: Found a bug? Please report it via [GitHub Issues](https://github.com/ScilifelabDataCentre/divbase/issues)
- **Security vulnerabilities**: Please report any security issues privately following the instructions in our [Security Policy](SECURITY.md)
- **Documentation**: Help improve our documentation by fixing typos, clarifying instructions, or adding new guides
- **Feature Requests**: Have an idea for a new feature? Consider opening an issue first so we can discuss it.

## Process

### 1. Fork and clone the repository

Start by [forking the DivBase repository](https://github.com/ScilifelabDataCentre/divbase/fork) to your GitHub account.

Clone your fork and create a branch for your changes:

```bash
git clone git@github.com:YOUR_USERNAME/divbase.git
cd divbase
git checkout -b your-feature-branch
```

### 2. Set up your development environment

Follow our [Developer Setup Guide](development/setup.md) to get your local development environment running:

### 3. Make your changes

- Write your code following our coding standards
- Add tests for new functionality
- Update documentation if needed
- Ensure all tests pass: `pytest -s`

### 4. Commit and push

Write clear, descriptive commit messages:

```bash
git commit -m "brief description of what you added"
git push origin your-feature-branch
```

## Opening a pull request

1. Navigate to the [original DivBase repository](https://github.com/ScilifelabDataCentre/divbase)
2. Click "New Pull Request"
3. Select your branch and create the pull request
4. Check the option to "Allow edits from maintainers" for easier collaboration

### Pull request guidelines

**In your PR description, please include:**

- A clear description of what your changes do
- Why the changes are needed
- Any breaking changes
- Any relevant advice on how to test out the changes.

## Testing

- Write tests for new features and bug fixes
- Ensure all existing tests still pass
- Tests are located in the `tests/` directory
- Run tests with: `pytest -s`

## Documentation

- Update relevant documentation for new features
- User guides are in `docs/user-guides/`
- Developer documentation is in `docs/development/`
