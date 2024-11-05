#!/usr/bin/env -S emacs -Q --script

;;==========================================================
;; Configure org
;;==========================================================

(require 'org)
(require 'ob-python)
(require 'python)
(require 'project)

(org-babel-do-load-languages
 'org-babel-load-languages
 '((python . t)))

(setq org-confirm-babel-evaluate        nil
      python-indent-guess-indent-offset nil)

(let ((project-root-dir (project-root (project-current))))
  (setq python-shell-extra-pythonpaths
        (append
         (mapcar 'expand-file-name
                 (list
                  (concat project-root-dir "src")
                  (concat project-root-dir "test")))
         python-shell-extra-pythonpaths)))

;;==========================================================
;; Now load pyretechnics.org and evaluate its src blocks
;;==========================================================

(message "Evaluating all source code blocks in org/pyretechnics.org")

(let ((script-directory (file-name-directory load-file-name)))
  (cd script-directory)
  (find-file "pyretechnics.org")
  (org-babel-execute-buffer)
  (save-buffer))

;;==========================================================
;; Obligatory calling shell protection
;;==========================================================

(setq argv nil)
