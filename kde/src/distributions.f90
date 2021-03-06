!------------------------------------------------------------------------------
! MODULE: distributions
! 
!> @author
!> Pavel Junker
! 
! DESCRIPTION: 
!> a small collection of benchmark pdfs, all are technically 1d
! 
! REVISION HISTORY:
! 13.03.2021 - Initial Version
!------------------------------------------------------------------------------
module distributions
	implicit none

	integer, parameter, private :: DP = kind(0.D0)

	abstract interface
		real(kind=DP) pure function rn_to_r(x)
			import :: DP
			real(kind=DP), intent(in) :: x(:)
		end function rn_to_r
	end interface

contains
! ========================== normal distribution pdf ==========================

!                                       #                                       
!                                     #####                                     
!                                   #########                                   
!                                  ###########                                  
!                                  ###########                                  
!                                 #############                                 
!                                ###############                                
!                               #################                               
!                               #################                               
!                              ###################                              
!                             #####################                             
!                             #####################                             
!                            #######################                            
!                            #######################                            
!                           #########################                           
!                          ###########################                          
!                          ###########################                          
!                         #############################                         
!                        ###############################                        
!                        ###############################                        
!                       #################################                       
!                      ###################################                      
!                     #####################################                     
!                    #######################################                    
!                    #######################################                    
!                  ###########################################                  
!                 #############################################                 
!                ###############################################                
!              ###################################################              
!          ########################################################### 

	pure function gauss(x) result(out)
		real(kind=DP), intent(in) :: x(:)
		real(kind=DP) :: out
        	out = exp(-0.5 * norm2(x)**2)
	end function
	
! ========================== flat-ish top bell curve ==========================

!                                       #                                       
!                                   #########                                   
!                                  ###########                                  
!                                 #############                                 
!                                 #############                                 
!                                ###############                                
!                                ###############                                
!                               #################                               
!                               #################                               
!                              ###################                              
!                              ###################                              
!                              ###################                              
!                             #####################                             
!                             #####################                             
!                             #####################                             
!                            #######################                            
!                            #######################                            
!                            #######################                            
!                           #########################                           
!                           #########################                           
!                          ###########################                          
!                          ###########################                          
!                         #############################                         
!                        ###############################                        
!                        ###############################                        
!                       #################################                       
!                      ###################################                      
!                    #######################################                    
!                  ###########################################                  
!              ###################################################  

	pure function plateau(x) result(out)
		real(kind=DP), intent(in) :: x(:)
		real(kind=DP) :: out
		out = 1/(1+norm2(x)**4)
	end function
	
! ================================== bimodal ==================================

!                      #                                 #                      
!                     ###                               ###                     
!                    #####                             #####                    
!                    ######                           ######                    
!                   #######                           #######                   
!                   ########                         ########                   
!                  #########                         #########                  
!                  #########                         #########                  
!                  ##########                       ##########                  
!                 ###########                       ###########                 
!                 ############                     ############                 
!                 ############                     ############                 
!                #############                     #############                
!                ##############                   ##############                
!                ##############                   ##############                
!               ###############                   ###############               
!               ################                 ################               
!              #################                 #################              
!              ##################               ##################              
!              ##################               ##################              
!             ###################               ###################             
!             ####################             ####################             
!            #####################             #####################            
!            ######################           ######################            
!           #######################           #######################           
!           ########################         ########################           
!          ##########################       ##########################          
!         ############################     ############################         
!        ###############################################################        
!      ###################################################################  

	pure function bimodal(x) result(out)
		real(kind=DP), intent(in) :: x(:)
		real(kind=DP) :: out
		out = exp(-1.5*(norm2(x)-1.5)**2) + exp(-1.5*(norm2(x)+1.5)**2)
	end function
	
! ========================== comb-like (more nodes) ===========================

!                                       #                                       
!                                       #                                       
!                                      ###                                      
!                                      ###                                      
!                                      ###                                      
!                                      ###                                      
!                                      ###                                      
!                                     #####                                     
!                                     #####                                     
!                                     #####                                     
!                                     #####                                     
!                             #       #####       #                             
!                            ##       #####       ##                            
!                            ###     #######     ###                            
!                            ###     #######     ###                            
!                           ####     #######     ####                           
!                           #####    #######    #####                           
!                           #####    #######    #####                           
!                           #####    #######    #####                           
!                          ######   #########   ######                          
!                          #######  #########  #######                          
!                          #######  #########  #######                          
!                          ###########################                          
!                         #############################                         
!                         #############################                         
!                   #     #############################     #                   
!                 ####   ###############################   ####                 
!                ###############################################                
!               #################################################               
!              ###################################################  

	pure function comb(x) result(out)
		real(kind=DP), intent(in) :: x(:)
		real(kind=DP) :: out
		out = (0.65d0+0.35d0*cos(6.5d0*norm2(x)))*exp(-0.5d0*norm2(x)**2)
	end function
	
! =========================== comb-like (less nodes) ==========================

!                                       #                                       
!                                      ###                                      
!                                      ###                                      
!                                     #####                                     
!                                     #####                                     
!                                     #####                                     
!                                    #######                                    
!                                    #######                                    
!                                    #######                                    
!                                   #########                                   
!                                   #########                                   
!                                   #########                                   
!                                  ###########                                  
!                                  ###########                                  
!                                 #############                                 
!                           ###   #############   ###                           
!                          ###########################                          
!                         #############################                         
!                         #############################                         
!                        ###############################                        
!                        ###############################                        
!                       #################################                       
!                       #################################                       
!                      ###################################                      
!                     #####################################                     
!                     #####################################                     
!                    #######################################                    
!                  ###########################################                  
!              ###################################################              
!           ######################################################### 

	pure function comb2(x) result(out)
		real(kind=DP), intent(in) :: x(:)
		real(kind=DP) :: out
		out = (0.8d0+0.2d0*cos(5.0d0*norm2(x)))*exp(-0.5d0*norm2(x)**2)
	end function

! ============================== skewed (1d only) =============================

!                                        #                                      
!                                        #                                      
!                                        ##                                     
!                                        ##                                     
!                                        ##                                     
!                                        ##                                     
!                                        ###                                    
!                                        ###                                    
!                                        ###                                    
!                                        ###                                    
!                                       #####                                   
!                                       #####                                   
!                                       #####                                   
!                                       #####                                   
!                                       ######                                  
!                                       ######                                  
!                                       ######                                  
!                                       #######                                 
!                                       #######                                 
!                                       #######                                 
!                                       ########                                
!                                       ########                                
!                                       #########                               
!                                       #########                               
!                                      ###########                              
!                                      ############                             
!                                      ############                             
!                                      ##############                           
!                                      ###############                          
!                                     ###################

	pure function skewed(x) result(out)
		real(kind=DP), intent(in) :: x(:)
		real(kind=DP) :: out
		out = 1 / (1 + exp(-20*x(1))) / (1 + exp(3*x(1)))
	end function

! ============================= a bit wider skewed ============================

!                    #                                                          
!                   ###                                                         
!                   ####                                                        
!                   ####                                                        
!                  ######                                                       
!                  #######                                                      
!                  #######                                                      
!                  ########                                                     
!                  #########                                                    
!                  ##########                                                   
!                 ###########                                                   
!                 ############                                                  
!                 #############                                                 
!                 ##############                                                
!                 ###############                                               
!                 ###############                                               
!                 ################                                              
!                ##################                                             
!                ####################                                           
!                #####################                                          
!                ######################                                         
!                #######################                                        
!                #########################                                      
!               ###########################                                     
!               #############################                                   
!               ###############################                                 
!               ##################################                              
!              ######################################                           
!              ###########################################                      
!             ####################################################    

	pure function skewed2(x) result(out)
		real(kind=DP), intent(in) :: x(:)
		real(kind=DP) :: out
		out = 1 / (1 + exp(-10*(x(1)+2))) / (1 + exp(x(1)+2))
	end function

end module distributions
