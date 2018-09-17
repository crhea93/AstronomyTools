function show_image() {

    var arc_per_pix = document.getElementById('Telescope').value;

    var file = document.getElementById('file').files[0];
    var reader  = new FileReader();
    // it's onload event and you forgot (parameters)
    reader.onload = function(e)  {
        var image = document.createElement("img");
        // the result image data
        image.src = e.target.result;
        var canvas = document.createElement("canvas");
        document.body.appendChild(canvas);

        /*canvas.width  = image.src .width;
        canvas.height = image.src .height;*/

        /*var context = canvas.getContext("2d");

        context.drawImage(image.src, 0, 0);*/
        /*document.body.appendChild(image);*/
    };
    // you have to declare the file loading
    reader.readAsDataURL(file);
}
