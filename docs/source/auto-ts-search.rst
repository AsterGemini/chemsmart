###########################
 Automated TS check/queue
###########################

This MVP follows the ts-scout Loop.md pass-2 spec (critic-accepted PASS 2
FINAL). It adds **no new job types** and does **not** run MLIP/UMA.

Two top-level commands:

-  ``chemsmart check LOG --spec case.yaml`` — parse a Gaussian log and print
   the mode table plus hard rejects.
-  ``chemsmart queue ts|irc|opt`` — write a shell file of existing
   ``chemsmart sub gaussian …`` commands. Never submitted.

Scientific acceptance of a TS stays with the user (final OK). Every submit
stays with the user.

*************************
 check
*************************

.. code:: bash

   chemsmart check ts.log --spec case1.yaml
   chemsmart check ts.log --spec case1.yaml --force
   chemsmart check ircf.log --spec case3.yaml --ts-file ts.log

Prints nimag, imaginary frequency, charge/mult, the top-8 heavy-atom pair
projections (no cutoff, with distances), spec-bond ranks, and Zmax
(information only).

Hard rejects
============

-  ``nimag≠1`` — hard reject. ``--force`` does not override.
-  Primary spec bond not #1 among heavy-atom pairs with r ≤ 3.2 Å — reject
   that ``--force`` overrides. The rank is always printed.

Everything else is labelled ``CANDIDATE: judge mode``.

Route and charge/mult
=====================

Expected (charge/mult 0/1):

-  TS: ``# opt=(ts,calcfc,noeigentest,maxstep=5) freq MN15 def2svp``
-  IRC: ``# MN15 def2svp irc(calcfc,recalc=6,forward|reverse,maxpoints=512,maxcycle=128)``
-  opt: ``# opt freq MN15 def2svp``

Mismatches are flagged (the old Negative 1 route is a CANDIDATE with a
ROUTE FLAG). This command does not edit ``calc/chemsmart/zn5.yaml``.

IRC logs
========

-  Fewer than 5 IRC points: print ``suggest qrc ±`` and the
   ``chemsmart sub gaussian … qrc -a 0.5`` command. Suggestion only; not
   written to a queue file and not submitted.
-  IRC reached minima but the end geometry fails the YAML endpoint bond
   criteria: print ``connectivity mismatch: judge`` and suggest no QRC.

Endpoint-opt logs: require nimag=0 and check the YAML endpoint bonds.

*************************
 queue
*************************

.. code:: bash

   chemsmart queue -p zn5 -c 0 -m 1 -o ts.sh ts -f guess.xyz
   chemsmart queue -p zn5 -o irc.sh irc -f ts.log
   chemsmart queue -p zn5 -o opt.sh opt --ircf ircf.log --ircr ircr.log

Writes a shell file. **Does not submit.** The TS line uses
``--additional-opt-options maxstep=5``. IRC is one ``irc`` command (forward
and reverse together). Endpoint opts are two ``opt`` lines in one file.

QRC is never queued here.

*************************
 Case spec YAML
*************************

The user writes the YAML (reaction-coordinate bonds and endpoint criteria).
Examples:

-  ``chemsmart/settings/templates/tscheck/case1.yaml`` — Case 1 TS1, C8–N29
-  ``chemsmart/settings/templates/tscheck/case1_ts2.yaml`` — Case 1 TS2, O1–C8
-  ``chemsmart/settings/templates/tscheck/case3.yaml`` — Case 3 TS2, C113–O114
   (INT2 < 1.5 Å, product > 2.3 Å)

*************************
 What stayed manual
*************************

-  Choosing the case / mechanism / ``--bond`` list in the YAML
-  Judging the printed mode table
-  Every ``chemsmart sub`` / running the queue file
-  Final OK on a verified TS
-  Cluster ``zn5.yaml`` (stale copy in ts-scout is a note for Chi, not edited)

UMA/MLIP, early-kill, and custom Gaussian templates are out of scope.
UMA is deferred until an HF token exists and a 3-TS mode-agreement test
passes.

Loop.md notes: ``docs/auto_ts_search_loop.md``.
