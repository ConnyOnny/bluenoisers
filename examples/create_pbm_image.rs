extern crate bluenoisers;
use bluenoisers::blue_noise;

fn main() {
    // create the samples
    let samples = blue_noise(vec![320., 240.], 4., 30);
    // create a boolean image to put the samples as pixels
    let mut image = vec![vec![false; 320]; 240];
    for s in samples {
        image[s[1] as usize][s[0] as usize] = true;
    }
    // write the image in portable bitmap format to stdout
    println!("P1");
    println!("{} {}", 320, 240);
    for y in 0..240 {
        for x in 0..320 {
            print!("{} ",
                   if image[y][x] {
                       "1"
                   } else {
                       "0"
                   });
        }
        println!("");
    }
}
