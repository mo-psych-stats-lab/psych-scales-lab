/* Firebase web configuration (optional).
   Leave as null to run the site without accounts; everything else keeps working and progress stays in the browser.
   To enable email/password accounts, paste the config object from:
   Firebase console > Project settings > General > Your apps > Web app > SDK setup and configuration > Config
   This object is NOT a secret (it identifies the project in every browser), but protect it as described in
   SETUP-FIREBASE.md: restrict the API key to your GitHub Pages domain and deploy firestore.rules.

   window.FIREBASE_CONFIG = {
     apiKey: "...",
     authDomain: "your-project.firebaseapp.com",
     projectId: "your-project",
     appId: "..."
   };
*/
window.FIREBASE_CONFIG = null;
/* Optional: App Check with reCAPTCHA v3 to block abusive scripts (see SETUP-FIREBASE.md, step 6). */
window.FIREBASE_APPCHECK_SITE_KEY = null;
