#!/usr/bin/env -S emacs -Q --script

;;==========================================================
;; Load all source code files within the repository
;; and detangle each one back into pyretechnics.org
;;==========================================================

(require 'org)

;; Prevent extra whitespace from being injected into code blocks
(setq org-src-preserve-indentation      t
      python-indent-guess-indent-offset nil)

(let* ((script-directory     (file-name-directory load-file-name))
       (root-directory       (file-name-parent-directory script-directory))
       (detangle-file-regexp (string-join '("\\.py$"
                                            "\\.pyx$"
                                            "\\.sh$"
                                            "\\.el$")
                                          "\\|"))
       (files-to-detangle    (seq-remove (lambda (file-name) (string-match-p "/org/\\|/doc/" file-name))
                                         (directory-files-recursively root-directory detangle-file-regexp))))
  ;; Enter the org/ directory and set the current buffer to pyretechnics.org so that (save-buffer) works later
  (cd script-directory)
  (find-file "pyretechnics.org")
  (mapc (lambda (source-file)
          (message "Detangling %s" source-file)
          (org-babel-detangle source-file)
          (save-buffer)) ; Save pyretechnics.org after each detangling operation
        files-to-detangle))

;;==========================================================
;; Obligatory calling shell protection
;;==========================================================

(setq argv nil)
