#!/usr/bin/env -S emacs -Q --script

;;==========================================================
;; Load Pyretechnics.org and tangle its code blocks to disk
;;==========================================================

(require 'org)
(require 'python-mode)

(find-file "Pyretechnics.org")
(org-babel-tangle)

;;==========================================================
;; Obligatory calling shell protection
;;==========================================================

(setq argv nil)
