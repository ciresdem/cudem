# CUDEM Developer Environment Setup & Contribution Workflow

This document describes how to set up a local CUDEM development environment, connect a GitHub fork, keep your code in sync with the main project, and submit pull requests (PRs) to the CUDEM repository.

It is intended for developers working on CUDEM modules (e.g., `blend.py`) or contributing new tools, enhancements, or documentation.

---

## 1 Configure Git (one time per machine)

Set your Git identity so commits appear correctly on GitHub:

```bash
git config --global user.name "Your Name"
git config --global user.email "your_email@example.com"
```

Use the **same email** associated with your GitHub account.

---

## 2 Fork the CUDEM Repository (do once)

1. Visit the official CUDEM repository:  
   **https://github.com/ciresdem/cudem**
2. Ensure you are signed into GitHub.
3. Click **Fork** (top-right).

Your fork will be created at:

```
https://github.com/<your-username>/cudem
```

---

## 3 Clone *Your Fork* Locally

```bash
git clone https://github.com/<your-username>/cudem.git
cd cudem
```

This creates a local working copy linked to your fork.

---

## 4 Add the Upstream CUDEM Repository

To receive updates from the main CUDEM project, add an upstream remote:

```bash
git remote add upstream https://github.com/ciresdem/cudem.git
```

Verify setup:

```bash
git remote -v
```

Expected:

```
origin    https://github.com/<your-username>/cudem.git
upstream  https://github.com/ciresdem/cudem.git
```

- **origin** → your fork  
- **upstream** → official CUDEM repo  

---

## 5 Work on the `dev` Branch

CUDEM development occurs in the `dev` branch.

Switch to it:

```bash
git checkout dev
```

Always sync before making changes:

```bash
git fetch upstream
git merge upstream/dev
```

Optional but recommended:

```bash
git push origin dev
```

This keeps your fork’s `dev` branch updated.

---

## 6 Update or Replace Files (example: modifying blend.py)

Overwrite an existing file:

```bash
cp /path/to/your/blend.py cudem/grits/blend.py
```

Check what changed:

```bash
git status
```

---

## 7 Commit and Push Your Changes

Stage the modified file:

```bash
git add cudem/grits/blend.py
```

Commit with a clear message:

```bash
git commit -m "Improve blend.py with slope-gated blending and parameter updates"
```

Push to your fork:

```bash
git push origin dev
```

---

## 8 Open a Pull Request (PR)

In your browser:

1. Go to your fork:  
   **https://github.com/<your-username>/cudem**
2. GitHub will suggest:  
   **"Compare & pull request"**

Set PR parameters:

| Setting             | Value                |
|---------------------|----------------------|
| **Base repository** | `ciresdem/cudem`     |
| **Base branch**     | `dev`                |
| **Head repository** | `<your-username>/cudem` |
| **Head branch**     | `dev`                |

Write a brief description of the changes and click **Create Pull Request**.

---

## 9 Keeping Your Fork Updated Over Time

To update your fork with changes from the main repository:

```bash
git checkout dev
git fetch upstream
git merge upstream/dev
git push origin dev
```

This prevents merge conflicts and ensures your development environment matches current CUDEM code.

---

## 10 Developer Install (Editable Mode)

To run CUDEM directly from your local source tree:

```bash
pip install -e .
```

This installs console entry points (e.g., `grits`) and ensures that edits in the source directory are immediately reflected in your environment.

---

## ✔ Summary of the Contribution Workflow

1. Fork → Clone → Add upstream  
2. Switch to `dev`  
3. Sync with upstream  
4. Make code changes  
5. Commit and push to your fork  
6. Open PR against `ciresdem/cudem:dev`  

This standard open-source workflow ensures clean version control and smooth collaboration across contributors.

---

## Suggestions and Extensions

Contributors are encouraged to enhance this document with:

- Environment setup instructions (conda or virtualenv)  
- Testing workflows  
- Style guidelines  
- Example CI-relevant commands  
- Module-specific development notes  

Additional PRs improving documentation are welcome.
