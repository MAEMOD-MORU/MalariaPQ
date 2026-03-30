$(document).ready(function() {
  // Hide the button initially
  $('#goTopButton').hide();

  // Show the button when the user scrolls past 500px
  $(window).scroll(function() {
    if ($(this).scrollTop() > 400) {
      $('#goTopButton').fadeIn();
    } else {
      $('#goTopButton').fadeOut();
    }
  });

  // Scroll to top when the button is clicked
  $('#goTopButton').click(function() {
    $('html, body').animate({ scrollTop: 0 }, 'slow');
    return false;
  });
});
