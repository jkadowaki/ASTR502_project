SELECT
  p.ra, p.dec, p.u, p.g, p.r, p.i, p.z
FROM PhotoObjAll AS p
WHERE
  p.ra BETWEEN 161.3375 AND 163.3375
  AND p.dec BETWEEN 50.05 AND 52.05
AND dbo.fDistanceArcMinEq(162.3375, 51.05, p.ra, p.dec) < 30.
