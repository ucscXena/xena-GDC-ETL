Contributing Guide
================================================================================

How you can contribute?
--------------------------------------------------------------------------------

- If you have encountered any bug while using this package, you can help us by
  submitting an issue in the GitHub issue tracker.
- If you have any feature in mind that is missing, you can also submit that in
  the issue tracker.
- You can also submit pull requests for the existing issues. But, please make 
  sure there is no open PR for the same issue.

How to contribute?
--------------------------------------------------------------------------------

- First fork the repository in your namespace.
- Now clone the repository locally by running:  

  ``git clone https://github.com/[your_username]/xena-GDC-ETL``

- Now create a new branch locally:

  ``git checkout -b new_branch``

- Make the required changes and stage them:

  ``git add .`` or ``git add path/to/files``

- Commit the staged changes:

  ``git commit -m "Some suitable commit message"``

  If there is an issue that your commit solves, also include that in the commit
  message, i.e. ``Closes #issue_no``. This will close the issue once your commit
  is merged.

- Now push the changes into remote:

  ``git push origin new_branch``

- And finally submit a Pull Request.

.. note::

  - Follow `PEP8 <https://pep8.org>`_ guidelines while writing Python codes.
    You can check if your code follows PEP8 guidelines by running ``flake8``.
  - Add docstrings for classes, functions and modules. We follow `Google's
    docstring guide <https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings>`_
    for them.
  - If you are adding tests, make sure to mark it with CI metadata if the test
    requires an internet connection. Read more
    `here <https://docs.pytest.org/en/latest/example/markers.html>`_.
      - Now if you want to execute tests that require internet connection, you can
        run them by:

      ``pytest -v -m CI``

      - And if you want to execute tests when you are offline, you can run them by:

      ``pytest -v -m "not CI"``

      - You can run all the tests simply by:

      ``pytest``
