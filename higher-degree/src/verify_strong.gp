\\ verify_strong.gp -- check, in exact arithmetic, every exclusion made by the STRONG
\\ tests of enum6.c (lines "NS 0 a1..ad", "RIGM n a1..ad", "RIGP n a1..ad").
\\ For irreducible f with V = Var(x_1):
\\   NS   : the exact minimum of Var on Z[x_1] \ Z must be < V;
\\   RIGM : alpha = x - n must have tau < 3/2 V and sqshape(f, alpha, V) = [];
\\   RIGP : alpha = n - x likewise.
\\ Prints any disagreement; returns counts [checked, reducible skipped, disagreements].
vstrong(infile, part, nparts) =
{
  my(L = externstr(concat("cat ", infile)), chk = 0, red = 0, bad = 0);
  for (i = 1, #L,
    if ((i - 1) % nparts != part, next);
    my(w = strsplit(L[i], " "), tag = w[1], n = eval(w[2]), c = apply(eval, w[3..#w]), f, d, V, ok);
    if (tag == "OUT", next);
    f = Pol(concat([1], c)); d = poldegree(f);
    if (!polisirreducible(f), red++; next);
    V = ((d-1)*c[1]^2/d - 2*c[2]) / d;
    if (tag == "NS", ok = (onu(ordZx(f)) < V),
      my(al = if (tag == "RIGM", x - n, n - x), ta = trace(Mod(al, f)) / d, s);
      s = sqshape(f, al, V);
      ok = (ta < 3/2*V) && !(s === -1) && (#s == 0));
    chk++;
    if (!ok, bad++; print("DISAGREE: ", L[i])));
  [chk, red, bad];
}
