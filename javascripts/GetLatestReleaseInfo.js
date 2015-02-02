function GetLatestReleaseInfo() {
	$.getJSON("https://api.github.com/repos/naturalis/supersmart/tags").done(
			function(json) {

				var release = json[0];
				var name = release.name;

				var zipURL = "https://github.com/naturalis/supersmart/archive/"
						.concat(name, ".zip");
				var tarURL = "https://github.com/naturalis/supersmart/archive/"
						.concat(name, ".tar.gz");
				$("#zipdownload").attr("href", zipURL);
				$("#tardownload").attr("href", tarURL);
				var vnumber = name.slice(1, name.length);
				$('#dirname').html("supersmart-".concat(vnumber));
			});
}
