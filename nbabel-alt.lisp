;; This version uses general-purpose vector functions based on simple
;; lists, and decoupled integration algorithms to allow easy
;; comparison of integration schemes.

(declaim (optimize (speed 3) (debug 0) (safety 0))
         (inline v+v v-v vlen^2 vlen s*v force))

(setf *read-default-float-format* 'double-float)

(defconstant +dt+ 1d-3)

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

(defgeneric load-stdin (c))

(defgeneric timestep (c integrator))

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

(defun stdin->strlist ()
  "read stdin as a list of strings"
  (loop for line = (read-line *standard-input* nil)
     while line
     collect line))

(defun split-spaces (s &optional acc)
  "tokenize a string by splitting at space characters"
  (let ((o (search " " s)))
    (if (null o)
        (reverse (remove-if (lambda (x) (string= "" x)) (cons s acc)))
        (split-spaces (string-trim '(#\ ) (subseq s o)) (cons (subseq s 0 o) acc)))))

(defmethod load-stdin ((c cluster))
  "initialize the cluster from stdin"
  (setf (stars c)
        (loop for line in (remove-if (lambda (x) (< (length x) 3))
                                     (mapcar #'split-spaces (stdin->strlist)))
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

(defmethod timestep ((c cluster) integrator)
  "update `CLUSTER` for one timestep using `INTEGRATOR` as the integration algorithm"
  (loop for i below (length (stars c))
     as star = (elt (stars c) i)
     do (multiple-value-bind (nx nv) (funcall integrator i)
          (setf (pos star) nx)
          (setf (vel star) nv))))

(defun force-on-particle (cluster)
  "create a closure around `CLUSTER` for calculating the force on a particle"
  (lambda (i p)
    "calculate the total force acting on particle `I`, at position `P`"
    (loop for star in (stars cluster)
       as xi = (pos star)
       as j = 0 then (incf j)
       as mj = (mass (elt (stars cluster) j))
       as fi = (if (not (= i j))
                   (force p xi mj)
                   (list 0 0 0))
       as acc = fi then (v+v acc fi)
       finally (return acc))))

(defun euler (c)
  "create a lambda representing the explicit euler algorithm. `C` is the star cluster."
  (let ((f (force-on-particle c)))
    (lambda (i)
      (let* ((star (elt (stars c) i))
             (xi (pos star))
             (vi (vel star))
             (xi+1 (v+v xi (s*v +dt+ vi))))
        (values
         xi+1
         (v+v vi (s*v +dt+ (funcall f i xi))))))))

(defun modified-euler (c)
  "create a lambda representing the modified euler algorithm. `C` is the star cluster."
  (let ((f (force-on-particle c)))
    (lambda (i)
      (let* ((star (elt (stars c) i))
             (xi (pos star))
             (vi (vel star))
             (vi+1 (v+v vi (s*v +dt+ (funcall f i xi))))
             (xi+1 (v+v xi (s*v +dt+ vi+1))))
        (values xi+1 vi+1)))))

(defun midpoint-euler (c)
  "create a lambda representing the midpoint euler algorithm. `C` is the star cluster."
  (let ((f (force-on-particle c)))
    (lambda (i)
      (let* ((star (elt (stars c) i))
             (xi (pos star))
             (vi (vel star))
             (a (funcall f i xi))
             (vi+1/2 (v+v vi (s*v (* 1/2 +dt+) a)))
             (vi+1 (v+v vi (s*v +dt+ a)))
             (xi+1 (v+v xi (s*v +dt+ vi+1/2))))
        (values xi+1 vi+1)))))

(defun leapfrog (c)
  "create a lambda representing the leapfrog algorithm. `C` is the star cluster."
  (let ((f (force-on-particle c)))
    (lambda (i)
      (let* ((star (elt (stars c) i))
             (xi (pos star))
             (vi (vel star))
             (ai (funcall f i xi))
             (xi+1 (v+v xi (s*v +dt+ vi) (s*v (* 1/2 +dt+ +dt+) ai)))
             (ai+1 (funcall f i xi+1)))
        (values
         xi+1
         (v+v vi (s*v (* 1/2 +dt+) (v+v ai ai+1))))))))

(defun main ()
  (let ((brk "=======================")
        (cluster (make-instance 'cluster))
        (tend 1d0))
    (load-stdin cluster)
    (multiple-value-bind (e0 pe ke) (energy cluster)
      (lineout "Time" "Kinetic Energy" "Potential Energy" "Total Energy" "Energy Error")
      (lineout brk brk brk brk brk)
      (summary 0d0 e0 ke pe e0)
      (let ((integrator (leapfrog cluster)))
        (loop for time from 0d0 to tend by +dt+
           as inc = 0 then (1+ inc)
           do (timestep cluster integrator)
             (when (zerop (mod inc 10))
               (multiple-value-bind (e pe ke) (energy cluster)
                 (summary (+ (* 10 +dt+) time) e ke pe e0))))))))

