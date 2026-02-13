/*  
Validates the user password on the registration form against a list of 10000 most common passwords.
*/

const passwordInput = document.getElementById('password');
const confirmPasswordInput = document.getElementById('confirm_password');
const registerForm = document.getElementById('registerForm');
let commonPasswords = new Set();

async function setupPasswordValidation() {
    try {
        const response = await fetch('/static/files/pwdb_top-10000.json');
        const data = await response.json();
        commonPasswords = new Set(data);

        passwordInput.addEventListener('input', validatePasswords);
        confirmPasswordInput.addEventListener('input', validatePasswords);

        registerForm.addEventListener('submit', (event) => {
            validatePasswords();
            if (!registerForm.checkValidity()) {
                event.preventDefault(); // Stop registration form submission
            }
        });

    } catch (err) {
        console.error("Failed to load list of common passwords", err);
    }
}

function validatePasswords() {
    const userPassword = passwordInput.value;
    if (commonPasswords.has(userPassword.toLowerCase())) {
        passwordInput.setCustomValidity("That password is one of the 10,000 most common passwords. Please choose a better one.");
    } else {
        // setting validity to empty string means no error. 
        passwordInput.setCustomValidity("");
    }
    if (confirmPasswordInput.value !== userPassword) {
        confirmPasswordInput.setCustomValidity('Passwords do not match');
    } else {
        confirmPasswordInput.setCustomValidity('');
    }
}

setupPasswordValidation();