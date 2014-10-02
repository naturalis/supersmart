var app = {
	/* set up the UI behaviours */
	initialize : function() {
      
		// add the handlers that scroll to the wanted section
		$('.scrollable').each(function(){
			var href = this.hash;
			$(this).on('click',function(){
				$('html, body').animate({
					scrollTop: $(href).offset().top
				}, 500);
				return false;
			});      
		});		
	}
}
