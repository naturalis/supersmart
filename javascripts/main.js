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
		
		// set the OS version for UI elements
		var OSName="Unknown OS";
		if (navigator.appVersion.indexOf("Win")!=-1) OSName="Windows";
		if (navigator.appVersion.indexOf("Mac")!=-1) OSName="MacOS";
		if (navigator.appVersion.indexOf("X11")!=-1) OSName="UNIX";
		if (navigator.appVersion.indexOf("Linux")!=-1) OSName="Linux";
		this.os = OSName;
		
		// set the terminal commands
		$('.term').each(function(){
			if ( OSName !== 'Windows' ) {
				$(this).text('$');
			}
		});		
	    
	    // load function to get link to latest SUPERSMART release				
		$(document).ready(function () {
			GetLatestReleaseInfo();  
		});  
		
		// set the toggle for shell comments
		$('.comment').each(function(){
			var comment = $(this);
			comment.slideUp(0);
			var button = $('<a title="Typical shell feedback">+</a>').bind('click',function() {
				if ( $(this).text() == '+' ) {
					comment.slideDown();
					$(this).text('-');
				}
				else {
					comment.slideUp();
					$(this).text('+');
				}
			});
			var widget = $('<span class="toggle"></span>');
			widget.append('[',button,']<br />');
			$(widget).insertBefore(comment);
		});					
	}
}
