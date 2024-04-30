#!/usr/bin/env -S emacs -Q --script

;;==========================================================
;; Load all files under the specified root-directory
;; and detangle each one back into pyretechnics.org
;;==========================================================

(require 'org)

;; Prevent extra whitespace from being injected into code blocks
(setq org-src-preserve-indentation t)

(defun valid-directory-p (directory-name)
  (and directory-name (file-directory-p directory-name)))

(defvar detangle-file-regexp (string-join '("\\.txt$"
                                            "\\.org$"
                                            "\\.json$"
                                            "\\.el$"
                                            "\\.sh$"
                                            "\\.py$"
                                            "\\.pyx$")
                                          "\\|"))

(let ((root-directory (car argv)))
  (if (not (valid-directory-p root-directory))
      (princ "Usage: detangle.el <source-directory>\n")
    (let ((script-directory  (file-name-directory load-file-name))
          (files-to-detangle (seq-remove (lambda (file-name) (string-match-p "/target/\\|/classes/" file-name))
                                         (directory-files-recursively root-directory detangle-file-regexp))))
      ;; Enter the org/ directory and set the current buffer to pyretechnics.org so that (save-buffer) works later
      (cd script-directory)
      (find-file "./pyretechnics.org")
      (mapc (lambda (source-file)
              (org-babel-detangle source-file)
              (save-buffer)) ; Save pyretechnics.org after each detangling operation
            files-to-detangle))))

;;==========================================================
;; Obligatory calling shell protection
;;==========================================================

(setq argv nil)
