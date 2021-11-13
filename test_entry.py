from io import StringIO
from pylint.reporters import text
from pylint.lint import Run


def main():
    pylint_opts = ['--reports=n', 'test_ternary.py']
    pylint_output = StringIO()
    reporter = text.TextReporter(pylint_output)
    Run(pylint_opts, reporter=reporter, do_exit=False)
    print(pylint_output.read()) # Show messages emitted by this pylint run


if __name__ == '__main__':
    main()
