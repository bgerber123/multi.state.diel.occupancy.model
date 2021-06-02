#multi-state likelhood, which is used for WAIC or CPO or k-fold cv
multi.state.likelihood=function(psi.matrix,det.matrix,data){
  
  pDay=det.matrix[,1]
  pNight=det.matrix[,2]
  pND.0=det.matrix[,3]
  pND.D=det.matrix[,4]
  pND.N=det.matrix[,5]
  pND.ND=det.matrix[,6]

  psi0=psi.matrix[,1]
  psiDay=psi.matrix[,2]
  psiNight=psi.matrix[,3]
  psiND=psi.matrix[,4]
  
#need to do this for each site
site.max.state=max(data,na.rm = TRUE)
tab.site=table(data)  

states.possible=c(0,0,0,0) #add frquency in each state, none, day, night, DN
names(states.possible)=c("1","2","3","4")

index=match(rownames(tab.site),names(states.possible))
states.possible[index]=as.integer(tab.site)


# Site could be not occupied or in any other state
if(site.max.state==1){ # 
  likelihood= psi0+     #not occupied
    psiDay*(1-pDay)^states.possible[1]+ #occupied in day state but not detected 
    psiNight*(1-pNight)^states.possible[1]+ #occupied in the night state but not detected 
    psiND*pND.0^states.possible[1] #occupied in state 4 but not detected
}
  
# Site could be in state 2 or in state 4. Not state 3
if(site.max.state==2){ 
 
  likelihood= psiDay*pDay^states.possible[2]*(1-pDay)^(states.possible[1])+ #in state 2 and detected and not detected
    psiND*pND.0^(states.possible[1])*pND.D^(states.possible[2])*pND.N^(states.possible[3])*pND.ND^(states.possible[4])    #in state 4 and not detected
    }

#Site could be in state 3 or in state 4
if(site.max.state==3 & states.possible[2]==0){ 

  likelihood= psiNight*pNight^states.possible[3]*(1-pNight)^(states.possible[1])+ #in state 3 and detected and not detected
    psiND*pND.0^(states.possible[1])*pND.D^(states.possible[2])*pND.N^(states.possible[3])*pND.ND^(states.possible[4])     #in state 4 and not detected
}

#Site has to be in state 4  
if(site.max.state==4 | (states.possible[2]>0 & states.possible[3]>0)){ 
  
  likelihood= psiND*pND.0^(states.possible[1])*pND.D^(states.possible[2])*pND.N^(states.possible[3])*pND.ND^(states.possible[4])
  }
likelihood #send this as output from the function
} #end function