# Entscheidungen vor einer PRL-Einreichung

## Kernaussage des aktuellen Entwurfs

Vlads Daten zeigen eine deutliche abstandsabhängige Reorganisation der winkelaufgelösten Zustandsdichte in stark komprimierten, geometrisch festgehaltenen Graphenstrukturen. Emils Modell illustriert einen möglichen geometrischen Hybridisierungsmechanismus. Ein topologischer Bulk-Gap, eine Chern-Zahl, echte spektrale Fraktalität oder ein quantisierter Pumpmechanismus sind damit noch nicht nachgewiesen.

Der Entwurf ist vollständig formuliert und gesetzt. Er ist bewusst als Author-review draft gekennzeichnet: Ein hochwertiger Satz ersetzt keine fehlenden Daten und Nachweise. Die derzeitige rein qualitative Aussage erscheint für PRL noch nicht stark genug. Der wichtigste nächste Schritt ist ein belastbarer, neuer physikalischer Nachweis, nicht eine schärfere Formulierung derselben Bilder.

## 1. Von Vlad benötigt: quantitative DFT-Grundlage

- CP2K-Ausgaben, Eigenwerte, PDOS/LDOS und Plot-Skripte für die gezeigten Panels; alle Energieumrechnungen und Normierungen.
- Bestätigung der Energieeinheit (die Bilder zeigen nur E - EF), Gaussian-Breite, Fermi-Referenz, Spin-Summation und des ausgegebenen unbesetzten Zustandsbereichs.
- Vollständige Inputs aller Abstände und der DZVP-Rechnungen; genaue Basisnamen und CP2K-Version.
- Behandlung der nicht konvergierten weißen Spalten. Diese nicht als Bandlücken auslegen und nicht interpolieren.
- Ein gezielter Konvergenztest mit größerer Basis, Grid-Cutoff, k-Punkten und Vakuumabstand.
- Mindestens ein Größen-/Randtest und eine Projektion auf das Innere der bedeckten unteren Schicht. Die obere C- und H-Zahl ändert sich deutlich mit dem Winkel.
- Geometrisch passende Monolayer-Referenzen mit gemeinsamem Energiebezug. Die obere H-terminierte Insel ist nicht dieselbe Referenz wie die untere periodische Schicht.
- Kräfte bzw. Relaxations-/Stabilitätsprüfung bei den sehr kleinen Abständen; keine Druckangabe allein aus dem Abstand ableiten.

## 2. Mit Emil klären: Modell und Topologie

- Die frühere Skala um 10^-14 stammt aus t_nn = exp(-1/0.03). Sie ist kein automatischer Nachweis eines numerisch unbrauchbaren Modells. Energie und DOS müssen aber gemeinsam umskaliert werden; dies ist nun im Supplement hergeleitet und numerisch getestet.
- Die relative Kopplung exp[(a-d)/xi] verändert sich beim Abstandsscan tatsächlich. xi zu ändern verändert zusätzlich die intralayer-Hopping-Reichweite.
- Festlegen, ob das Modell nur ein generischer quasiperiodischer Vergleich ist oder quantitativ an Graphen angepasst werden soll. Der aktuelle Entwurf verwendet nur die erste, belegbare Interpretation.
- Radius-/Broadening-Konvergenz und geeignete Bulk-/Randtrennung für ausgewählte Kandidatenlücken.
- Einen belastbaren phasonabhängigen Gap-Label-, Chern- oder Pumpnachweis für mindestens eine durchgehende Lücke erbringen. Ein DOS-Muster allein ist kein Topologienachweis; eine Winkelachse ist weder k-Achse noch Sliding-Zyklus.
- Bei einer IDS-Gap-Label-Auswertung prüfen, welche Koeffizienten entlang der gewählten Geometriefamilie tatsächlich getrennt identifizierbar sind. Insbesondere liefert ein konstanter Determinantenterm nicht automatisch eine eindeutige zweite Chern-Zahl.

## 3. Neuheit gegenüber vorhandener Literatur

Ni et al., Communications Physics 2, 55 (2019), DOI 10.1038/s42005-019-0151-7, ist die wichtige experimentelle Motivation: akustisches Butterfly-Spektrum plus phasonabhängige Randzustände. Rosa, Ruzzene und Prodan, Communications Physics 4, 130 (2021), DOI 10.1038/s42005-021-00630-3, ist die unmittelbarere theoretische Vorarbeit für verdrillte Bilayer.

Abstand bzw. Druck als Kontrollparameter für Graphen ist bereits bekannt, insbesondere Yu et al., Physical Review B 102, 045113 (2020), und Palamara et al., npj 2D Materials and Applications 9, 100 (2025), DOI 10.1038/s41699-025-00616-7. Die Neuheit muss daher präziser sein als nur „Kompression verändert das Spektrum“. Besonders überzeugend wäre die Verbindung eines atomistisch robusten elektronischen Gaps mit dem erwarteten höherdimensionalen topologischen Mechanismus.

## 4. Formales und Freigabe

- Autorenreihenfolge, Affiliationen, Corresponding Author, Förderhinweise und Beiträge bestätigen. Die momentane Reihenfolge Vladislav Efremkin, Emil Prodan, Thomas D. Kühne ist ein Vorschlag.
- Datenarchiv und Freigabe der gelieferten Daten vereinbaren; keine privaten E-Mail-Exporte veröffentlichen.
- AI-Hinweis nach Autorenprüfung und aktueller Verlagspolitik finalisieren.
- Alle derzeit im Manuskript offen benannten Reproduzierbarkeitsfragen vor einer Einreichung beantworten oder die Aussage entsprechend enger fassen.

Es wurden keine neuen DFT-Rechnungen durchgeführt und keine Nachrichten oder Einreichungen versandt. Neu gerechnet wurden ausschließlich die fünf kleinen Tight-Binding-Parameterreihen des gelieferten Modells.
