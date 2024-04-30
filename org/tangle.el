#!/usr/bin/env -S emacs -Q --script

;;==========================================================
;; Load pyretechnics.org and tangle its code blocks to disk
;;==========================================================

(require 'org)
(require 'python)

(setq python-indent-guess-indent-offset nil)

(let ((script-directory (file-name-directory load-file-name)))
  (cd script-directory)
  (find-file "pyretechnics.org")
  (org-babel-tangle))

;;==========================================================
;; Obligatory calling shell protection
;;==========================================================

(setq argv nil)
