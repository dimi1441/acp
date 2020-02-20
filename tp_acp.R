# Liste des taches à faire: Amélioration des graphiques.


# Pour la suite:
#
# - Mettre les titres des figures
# - Trouver une meilleure nuance de bleu pour les graphiques(comme celle du support qui nous sert d'inspiration)
# - Améliorer l'affichage des noms en bas des graphiques.
# - Dessiner les cercles de corrélation et inserer les flèches.
# - Colorier les points présents sur les cercles et les différents plans factoriels.


verification = function(x, d){
  if((is.matrix(x) | is.data.frame(x)) & (is.matrix(d) | is.data.frame(d))){ # ON sassure que les 2 sont desmatrices ou des dataframe
    if(nrow(x) == nrow(d)){ # On vérifie que les 2 possèdent le meme nombre de lignes
      if(a){ # On vérifie que les données sont bien numériques
        if(identical(is.na(x), matrix(data=FALSE, nrow(x), ncol(x))) & identical(is.na(d), matrix(data=FALSE, nrow(d), ncol(d)))){ # On verifie l'absence de valeurs manquantes.
          print("4")
          return (TRUE) 
        }
      }
    }
  }
  return (FALSE)
}

inverse = function(x){
  if(is.vector(x)){
    result = c(1)
    for(i in 1:length(x)){
      if(x[i] == 0){
        result[i] = 0
      }
      else{
        result[i] = 1/(x[i])
      }
    }
    return (result)
  }
}

carte_individus = function(donnees, x_label, y_label, noms, axes){
  plot(donnees[, 1], donnees[, 2],type ="p",xlab = x_label, ylab = y_label, main = paste("ACP des individus : plan factoriel", axes[1],"-", axes[2] ), pch=19)
  abline(h=0, v=0)
  text(donnees[, 1], donnees[, 2], labels = noms, cex = 0.65, pos=3)
}

cartes_individus = function(C, nombre_axes, noms){
  i = 1
  for(i in 1:nombre_axes){
    j = i+1
    while(nombre_axes >= j){
      carte_individus(cbind(C[, i], C[, j]), paste("Dim", i, sep = ""), paste("Dim", j, sep = ""), noms, c(i,j))
      j = j+1
    }
  }
}

plan_factoriel = function(donnees, x_label, y_label, noms, axes){
  a=seq(0,2*pi,length=100)
  plot(cos(a), sin(a), type='l', lty=3,xlab=x_label, ylab=y_label, main="Cercle des corrélations" )
  arrows(0, 0, donnees[, 1], donnees[, 2], col='blue')
  text(donnees[, 1], donnees[, 2], labels = noms, cex = 0.75)
}

plans_factoriels = function(C, nombre_axes, noms){
  i = 1
  for(i in 1:nombre_axes){
    j = i+1
    while(nombre_axes >= j){
      plan_factoriel(cbind(C[, i], C[, j]), paste("Dim", i, sep = ""), paste("Dim", j, sep = ""), noms, c(i,j))
      j = j+1
    }
  }
}

centrer_reduire = function(x, d){
  
  g = cbind(t(x)%*%d%*%cbind(rep(1, length=nrow(x))))
  y = round(x - cbind(rep(1, length=nrow(x))) %*% t(g), 3) # Centrer
  
  print("Centrée")
  print(y)

  v = round(t(y)%*%d%*%y, 4)
  print("v")
  print(v)

  d_1_sur_s = diag(inverse(sqrt(diag(v))))
  
  z = y%*%d_1_sur_s

  return (z)
}

kaiser = function(valeurs_propres, nombre_variables){
  inertie_totale = sum(valeurs_propres)
  inertie_moyenne = inertie_totale / nombre_variables 
  nombre_axes = length(valeurs_propres[valeurs_propres > inertie_moyenne])
  print("les nombre d'axes à selectionner d'apres la methode Kaiser est:")
  print(nombre_axes)
  return (nombre_axes)
}

qualites_de_representations = function(C, nombre_axes, noms, choice){
  if(choice == "var"){
    nature = "variables" # tableau de taille ncol(C) de 1 sur racine de p
  }
  if(choice == "ind"){
    nature = "individus" # tableau des racines de lambda de taille ncol(C)
  }
  qrs = cbind(rep(0, nrow(C)))
  for(i in 1:nombre_axes){
    qr = inverse(apply(C*C, 1, sum)) * cbind(C[,i]*C[,i])
    qr_in_row = t(qr)
    barplot(qr_in_row, col = "blue", names.arg=noms, main = paste("ACP des", nature, ": Qualite de representation sur dim", i), las=2)
    qrs = cbind(qrs, qr)
  }
  qrs = qrs[,-1] # On supprime la première colonne.
  return (qrs)
}

contributions = function(C, Poids, values, noms, choice){
  if(choice == "var"){
    threshold = rep(1/sqrt(p), nombre_axes) # tableau de taille ncol(C) de 1 sur racine de p
    nature = "variables"
  }
  if(choice == "ind"){
    threshold = c(1, p) # tableau des racines de lambda de taille ncol(C)
    nature = "individus"
  }
  contribs = cbind(rep(0, nrow(C)))
  for(i in 1:nombre_axes){
    contrib = (1/values[i])*cbind(diag(Poids)*cbind(C[,i]*C[,i]))
    contrib_percent = t(contrib)
    barplot(contrib_percent, col = "blue", names.arg=noms, main = paste("ACP des", nature, ": Contribution sur dim", i), las=2)
    abline(h= threshold[i], col = "red", lwd=3, lty=2)
    contribs = cbind(contribs, contrib)
  }
   # On supprime la première colonne.
  return (contribs[,-1])
}

valeurs_vecteurs = function(v, m){
  propres = eigen(v%*%m)
  nb = length(propres$values[propres$values > 0.001]) # Le nombre de valeurs propres positives
  
  valeurs_propres = propres$values[c(1:nb)]
  print("Les valeurs")
  print(valeurs_propres)
  
  vecteurs_propres = propres$vectors[,c(1:nb)]
  print("les vecteurs")
  print(vecteurs_propres)
  
  resultat = list(values=valeurs_propres, vectors = vecteurs_propres)
  return (resultat)
}

eboulis = function(valeurs){
  barplot(as.vector(valeurs), col = "blue", lty=3, xlab="indices", ylab="valeurs", main="Ebouli des valeurs propres")
}

ACP_variable = function(c, d, m, data, values, noms){
  # determiner les axes factoriels
  factoriels = c %*% diag(inverse(sqrt(values)))
  print("axes factoriels")
  print(factoriels)
  
  # determiner les composantes principales
  composantes = round(t(data) %*% d %*% factoriels, 4)
  print("composantes principales")
  print(composantes)
  
  # Calculer les contributions des variables aux axes
  contrib = contributions(composantes, m, values, noms, "var")
  print("Les contributions sur les axes")
  print(contrib)
  
  # Calcuer les qualités de représentation des variabless
  qr = qualites_de_representations(composantes, nombre_axes, noms, "var")
  print("Les qualités de représentation sur les axes")
  print(qr)
  
  # Dessins
  plans_factoriels(composantes, nombre_axes, noms)
}

ACP = function(x, d){
  # Si les vérifications suivantes ne sont pas toutes respectées, on sort du programme:
  if(TRUE){
    if(nrow(d) == ncol(d) & sum(diag(d)) == 1 & identical(diag(diag(d)), d)){ # 2) d est un matrice diagonale dont la somme des poids donne 1
      
      # le nombre de variables
      p = ncol(x)
      
      ind_names = row.names(x)
      var_names = names(x)
      
      m = matrix(data = c(1), nrow = ncol(x), ncol = ncol(x)) # 1)
      v = matrix(data = c(1) , ncol(x), ncol(x))
      
      x = as.matrix(x)
      nombre_axes = 3
      
      if (!(identical(apply(x, 1, sum), rep(0, length=nrow(x))))){ # Teste si x n'est pas centré
        x = centrer_reduire(x, d)
        print("Matrice centreé et réduite")
        print(x)
        print("Matrice des correlations")
        v = t(x)%*%d%*%x
        print(v)
        print("Matrice des métriques")
        m = diag(rep(1, ncol(x))) # La matrice identité
        print(m)
      }
      else{
        print("Matrice des variances_covariances")
        v = t(x)%*%d%*%x # matrice des variance-covariance
        print("Matrice des métriques")
        m = diag(inverse(diag(variance))) # m recoit D1/s au carré
        print(m)
      }
      
      # Vecteurs et valeurs propres
      propres = valeurs_vecteurs(v, m)
      values = propres$values
      vectors = propres$vectors
      
      eboulis(values)
      
      print("Les composantes princiales")
      c = -1*x%*%m%*%vectors # Lorsque l'on met moins ici, les resultats sont les meme qu'avec le fonction pca
      print(c)
      
      # Ici, on crée la carte, qui peut etre un sev de dim 1 ou 2.
      carte = c[1]
      if(is.vector(c) == TRUE){
        # On fait la meme chose que hors du if sauf qu'ici, on ne dessine pas de cercle de corrélation
      }
      
      print("taille colonne")
      print(nrow(c))
      print("nombre de noms de individus")
      print(length(ind_names))
      
      # Calculer les contributions des individus aux axes
      contrib = contributions(c, d, values, ind_names, "ind")
      print("Les contributions sur les axes")
      print(contrib)
      
      # Calcuer les qualités de représentation des individus
      qr = qualites_de_representations(c, nombre_axes, ind_names, "ind")
      print("Les qualités de représentation sur les axes")
      print(qr)
      
      # Dessins
      cartes_individus(c, nombre_axes, ind_names)
      
      ACP_variable(c, d, m, x, values, var_names)
      
    }
    else{
      print("La matrice d doit etre diagonale avec sa trace egale à 1.")
    }
  }
}

library("FactoMineR")
library("factoextra")
library("ggplot2")


data = decathlon2
data.active = data[, -ncol(data)]
x = data.active
print(x)
d = diag(rep(1/nrow(x), nrow(x)))

ACP(x, d)