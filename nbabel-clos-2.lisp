(ql:quickload :cl-ppcre :silent t)
(use-package :cl-ppcre)

;; this version is CLOS with general-purpose vector functions based on simple lists,
;;   and a decoupled leapfrog algorithm

(declaim (optimize (speed 3) (debug 0) (safety 0))
         (inline v+v v-v vlen^2 vlen s*v force))

(setf *read-default-float-format* 'double-float)

(defconstant dt 1d-3)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; some vector utility functions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun vlen^2 (r &optional (acc 0))
  "length of a vector, squared"
  (cond ((null r) acc)
        ((numberp r) (+ acc (* r r)))
        (t (vlen^2 (cdr r) (+ acc (* (car r) (car r)))))))

(defun vlen (r)
  "length of a vector"
  (sqrt (vlen^2 r)))

(defun s*v (s v)
  "scalar times vector"
  (loop for i in v
     collect (* s i)))

(defun v-v (&rest r)
  "vector subtraction, first vector minus all remaining vectors passed in"
  (flet ((rec (p1 p2)
           (loop for a in p1
              for b in p2
              collect (- a b))))
    (reduce #'rec r)))

(defun v+v (&rest r)
  "vector addition, sum of all vectors passed in"
  (flet ((rec (v1 v2)
           (loop for x in v1
              for y in v2
              collect (+ x y))))
    (reduce #'rec r)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; a class representing a star (i.e. particle)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass star ()
  ((pos :accessor pos :initarg :pos :initform (list 0d0 0d0 0d0))
   (vel :accessor vel :initarg :vel :initform (list 0d0 0d0 0d0))
   (mass :accessor mass :initarg :mass :initform 0d0)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; a class representing a star cluster
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defclass cluster ()
  ((stars :accessor stars :initarg :stars :initform '())))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; generic methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defgeneric potential-energy (c))

(defgeneric kinetic-energy (c))

(defgeneric energy (c))

(defgeneric load-file (c fname))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; specialized methods
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defmethod potential-energy ((c cluster))
  "potential energy calculation"
  (let* ((stars (stars c))
         (n (length stars)))
    (loop for i below n sum
         (loop for j from (1+ i) below n
            as star-i = (elt stars i)
            as star-j = (elt stars j)
            sum (/ (* -1 (mass star-i) (mass star-j)) (vlen (v-v (pos star-j) (pos star-i))))))))

(defmethod kinetic-energy ((c cluster))
  "kinetic energy calculation"
  (loop for star in (stars c)
     sum (* 1/2 (mass star) (vlen^2 (vel star)))))

(defmethod energy ((c cluster))
  "return total energy, potential energy, and kinetic energy"
  (let ((p (potential-energy c))
        (k (kinetic-energy c)))
    (values (+ p k) p k)))

(defun force (ri rj mj)
  "calculate the force between two particles"
  (let* ((p (v-v rj ri))
         (c (* mj (expt (vlen^2 p) -1.5))))
    (s*v c p)))

(defun file->strlist (filename)
  "read a text file as a list of strings"
  (with-open-file (stream filename)
    (loop for line = (read-line stream nil)
       while line
       collect line)))

(defmethod load-file ((c cluster) fname)
  "initialize the cluster from a file"
  (setf (stars c)
        (loop for line in (remove-if (lambda (x) (< (length x) 3))
                                     (mapcar (lambda (x) (split " +" x)) (file->strlist fname)))
           as parsed = (mapcar (lambda (x) (with-input-from-string (l x) (read l))) line)
           as star = (make-instance 'star
                                    :mass (elt parsed 1)
                                    :pos (apply #'list (subseq parsed 2 5))
                                    :vel (apply #'list (subseq parsed 5)))
        collect star)))

(defmethod print-object ((obj star) out)
  (print-unreadable-object (obj out :type t)
    (format out "(~a ~a ~a) (~a ~a ~a) ~a"
            (first (pos obj))
            (second (pos obj))
            (third (pos obj))
            (first (vel obj))
            (second (vel obj))
            (third (vel obj))
            (mass obj))))

(defmethod print-object ((obj cluster) out)
  (print-unreadable-object (obj out :type t)
    (progn
      (format out "cluster: ~%")
      (loop for star in (stars obj) do (format out " ~S~%" star)))))

(defun lineout (time ke pe total err)
  (format t "~A~24T~A~48T~A~72T~A~96T~A~%" time ke pe total err))

(defun summary (time tot ke pe e0)
  (lineout time ke pe tot (/ (- tot e0) e0)))

(defun main ()
  (let ((brk "=======================")
        (cluster (make-instance 'cluster))
        (tend (with-input-from-string (l (or (second (uiop:command-line-arguments)) "1")) (read l))))
    (load-file cluster (or (car (uiop:command-line-arguments)) "input128"))
    (flet ((leapfrog (f)
             "create a closure representing the leapfrog algorithm. `F` is a function that calculates a particle's instantaneous acceleration."
             (lambda (i)
               (let* ((star (elt (stars cluster) i))
                      (xi (pos star))
                      (vi (vel star)))
                 (let* ((ai (funcall f i xi))
                        (xi+1 (v+v xi (s*v dt vi) (s*v (* 1/2 dt dt) ai)))
                        (ai+1 (funcall f i xi+1)))
                   (values
                    xi+1
                    (v+v vi (s*v (* 1/2 dt) (v+v ai ai+1))))))))
           (Fi (i p)
             "calculate the total force acting on particle i, at position p"
             (loop for star in (stars cluster)
                as xi = (pos star)
                as j = 0 then (incf j)
                as mj = (mass (elt (stars cluster) j))
                as fi = (if (not (= i j))
                            (force p xi mj)
                            (list 0 0 0))
                as acc = fi then (v+v acc fi)
                finally (return acc))))
      (multiple-value-bind (e0 pe ke) (energy cluster)
        (lineout "Time" "Kinetic Energy" "Potential Energy" "Total Energy" "Energy Error")
        (lineout brk brk brk brk brk)
        (summary 0d0 e0 ke pe e0)
        (let ((integrator (leapfrog #'Fi)))
          (loop for time from 0d0 to tend by dt
             as inc = 0 then (1+ inc)
             do (loop for i below (length (stars cluster))
                   as star = (elt (stars cluster) i)
                   do (multiple-value-bind (nx nv) (funcall integrator i)
                        (setf (pos star) nx)
                        (setf (vel star) nv)))
               (when (zerop (mod inc 10))
                 (multiple-value-bind (e pe ke) (energy cluster)
                   (summary time e ke pe e0)))))))))

