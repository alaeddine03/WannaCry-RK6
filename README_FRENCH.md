Le script définir des fonctions pour analyser la propagation et l'impact du rançongiciel WannaCry en utilisant des modèles épidémiologiques. 
La première partie du code comprend des fonctions pour les dérivées des modèles épidémiologiques SIR (Susceptible, Infecté, Rétabli) et 
SIS (Susceptible, Infecté, Susceptible), qui sont couramment utilisés en épidémiologie pour modéliser la propagation des maladies.

Définitions de fonctions pour les dérivées du modèle SIR :

-La fonction dSIRdt calcule les taux de changement (dSdt, dIdt, dRdt) pour les compartiments de la population susceptibles (S), infectés (I) et rétablis (R).
-Les paramètres comprennent beta (taux d'infection), gamma (taux de rétablissement) et N (population totale).

Définitions de fonctions pour les dérivées du modèle SIS :

-La fonction dSISdt calcule les taux de changement (dSdt, dIdt) pour les compartiments susceptibles (S) et infectés (I).
-Le modèle SIS diffère du modèle SIR car il n'inclut pas de compartiment rétabli ; les individus passent directement de l'état infecté à l'état susceptible.
-Les mêmes paramètres (beta, gamma, N) sont utilisés.

Fonction d'estimation d'erreur :

-Cette partie du code semble inclure une fonction pour estimer l'erreur entre deux étapes temporelles, probablement utilisée dans des simulations ou 
des solutions numériques de ces modèles.
